---
layout: default
title: Samtools - Workflows
---

# Workflows

* [FASTQ to BAM/CRAM processing](#fastq_to_bam)
* [WES Mapping to Variant Calls - Version 1.0](#mapping_to_variant)
* [Filtering of VCF files](#variant_filtering)
* [Using CRAM within Samtools](#mapping_to_cram)

## <a name="fastq_to_bam"></a>FASTQ to BAM / CRAM - Version 1.0

Sequencing instruments produce unaligned data, typically in FASTQ
format.  It is possible to store unaligned data in BAM or CRAM, and
indeed it may be preferable as it permits meta-data in the header and
per-record auxiliary tags, however in this workflow we consider the
end product to be a sorted aligned BAM or CRAM file so we cover that
at the end.

There are two primary ways of producing this:

* Alignment / mapping to a known reference
* De-novo assembly

#### Alignment / mapping to a known reference

There are many tools for mapping sequences (reads) against a reference
file.  We do not recommend any single tool over another, but for
illustrative purposes are using [Minimap2](https://github.com/lh3/minimap2).

There are processing steps which together make a pipeline:

* Map / align
* Fix mate-pair issues
* Mark duplicates (part 1: preparation)
* Sort to positional order
* Mark duplicates (part 2: marking)
* Convert to final file format

We initially break these down into individual stages, reading and
writing from disk, but it is worth considering using a UNIX pipeline
to link these together avoiding temporary files.

#### Mapping

A basic mapping tool may take a single or pair of FASTQ files and
align against either a pre-built index or directly against a fasta file
(as used below).  The most naive approach is simply to save this file
as SAM:

    minimap2 -t 8 -a -x sr C.Elegans.fa SRR065390_1.fastq SRR065390_2.fastq -o CE.sam

This outputs in SAM (-a), uses 8 threads (-t 8), with options for
paired end short read (-x sr).

This output file will be in the original input file order, hence the
read pairs will be collated next to each other.  This is important as
the next step requires name-collated data.
Note some aligners may shuffle the data a bit when multi-threading is
enabled, but in all cases the output will still be name-collated.

Note if you use `nohup` to run a command and redirect stdout using
e.g. `> CE.sam`, then you may mix stdout and stderr together.  This
can be prevented by also adding `2>err` to explicitly save stderr to a
separate file or preferably use an explicit output option, such as the
`-o CE.sam` above.

#### Fixing of mate-pair issues

The `samtools fixmate` tool corrects any flaws in read-pairing that
may have been introduced by the aligner.  Sadly a number of them have
subtle bugs and quirks, so this can be considered as a proof-reading
step.  It ensures the SAM FLAG, RNEXT, PNEXT and TLEN fields are
correct and internally consistent.

Note `samtools fixmate`, as with other sub-commands, can read SAM
directly without needing an explicit option.  You may find older
tutorials which use `-S` to indicate SAM or have an explicit SAM to
BAM conversion using `samtools view`.  This is not necessary.

    samtools fixmate -O bam,level=1 CE.sam fixmate.bam

The `-O bam,level=1` requests the fastest level of BAM compression
for the output file.  We could also use `level=0` (or use `-u`) for
uncompressed output.

#### Marking duplicates (part 1: preparation)

Marking duplicates requires some analysis performed on data in
read-name collated order, and some performed in genome-position
order.  Rather than use a markdup tool that internally does sorting,
we break it down into two stages so it can be slotted in to our
pipeline without requiring any *additional* sorting requirements.

See [Duplicate Marking](algorithms/duplicate.html) for further details.

The first stage is also implemented with `samtools fixmate`, so we can
amend the previous step with an additional flag:

    samtools fixmate -O bam,level=1 -m CE.sam fixmate.bam

This adds mate cigar (MC) and mate score tags (ms) which will be used
later by `samtools markdup` proper.

#### Sorting to positional order

The data can now be converted to genome chromosome and coordinate
order with:

    samtools sort -l 1 -@8 -o pos.srt.bam -T /tmp/example_prefix fixmate.bam

Sort is highly parallel so the `-@8` option here enables to use of 8
additional CPU threads.  It can also be sped up by providing it with
more memory, but note the memory option (`-m`) is per-thread.  The `-l
1` indicates level 1 compression again.  We could also specify `-O
bam,level=1` as used above.

#### Marking duplicates (part 2: marking)

The main core of marking duplicates may now be ran on the
position-sorted file, utilising the extra tags added during the
`fixmate` step.

    samtools markdup -O bam,level=1 pos.srt.bam markdup.bam

#### Conversion to final file format

At this point you can convert to a more highly compressed BAM or to
CRAM with `samtools view`

    samtools view -@8 markdup.bam -o final.bam

or

    samtools view -T C.Elegans.fa -@8 markdup.bam -o final.cram

Note if there is no other processing to do after markdup, the final
compression level and output format may be specified directly in that
command.  See below for an example.

#### Joining it all together

You may notice that the output from each stage is used solely as the
input to the next stage.  We'd want to tidy up these intermediates at
the end too.  While the above may allow you to specify the commands in
a workflow language (e.g. CWL) they are not the most efficient method.

Using Unix pipelines is a faster approach.  Furthermore the lack of
needing temporary files on disk (excepting any internal temporary
files output by `sort`) means we can use uncompressed BAM for maximum
speed.  All samtools commands accept `-` as a synonym for stdin and stdout.
Most commands have a `-u` option to request uncompressed BAM,
but all should also accept a more explicit `-O bam,level=0`.

An example pipeline would be:

    minimap2 -t 8 -a -x sr C.Elegans.fa SRR065390_[12].fastq  | \
    samtools fixmate -u -m - - | \
    samtools sort -u -@2 -T /tmp/example_prefix - | \
    samtools markdup -@8 --reference C.Elegans.fa - final.cram

Adjusting the number of threads used by each step may be organism and
system specific, but samtools will only use as many as it needs so
specifying too many threads may not be too detrimental.

When estimating the total number of concurrent threads to allocate,
consider that the sort step is a crunch point that separates the steps
before it from the step afterwards.  The mapping, fixmate and partial
sort to temporary file steps will operate in parallel.  Once complete,
the sort merge (from temporary files) and markdup steps will then run
in parallel.

The alignment is typically the bottleneck, so some pipelines may wish
to split the input into multiple jobs and distribute multiple mapping
jobs.  These could be combined together (either before or after
`fixmate`) using `samtools cat`.

Finally we would recommend enabling `pipefail` in your shell if
available.  In Bash this is:

    set -o pipefail

Normally the shell will only report the exit status of the final
command in the pipeline.  With pipefail, a failure in any part of the
pipe will make the entire pipeline fail.

#### Converting back to FASTQ

None of the steps above removed data provided the aligner included the
unmapped reads rather than discarding them.  (Most will do this
automatically.)

The original FASTQ can be produced, albeit perhaps not in the exact
same order, by sorting your position-sorted back to name-sorted order
and then running `samtools fastq`.  The aligners will not care if the
data is *precisely* the same order as produced by the sequencing
instrument.   Storing your FASTQ in aligned BAM or CRAM format can be
one very simple way of reducing long-term storage requirements.

    samtools sort -n -@8 final.cram | \
    samtools fastq - -1 dat_1.fq -2 dat_2.fq > /dev/null

### De-novo assembly

We don't wish to cover this in detail as it is a complex topic and
readers should follow the procedures recommended by the authors of the
assembly software.

However note that most sequence assemblers produce a consensus, in
FASTA or FASTQ format and not individual alignments of every sequence.
If you wish to get the BAM or CRAM file of data aligned against this
consensus, for purposes of curation or downstream analysis, then
simply follow the Mapping section above (with an additional step, if
required, of building an index on the assembly consensus prior to
mapping).

Note if using CRAM then care must be taken to ensure you do not lose
the consensus file, or to amend the CRAM options to either embed the
reference into the CRAM file or to use non-reference based
compression.  An example of embedding the reference is:

        samtools view -O CRAM,embed_ref in.sam -o out.cram


## <a name="mapping_to_variant"></a>WGS/WES Mapping to Variant Calls - Version 1.0
The standard workflow for working with DNA sequence data consists of three major steps:

* Mapping
* Improvement
* Variant Calling
	

### Mapping
For reads from 70bp up to a few megabases we recommend using [BWA MEM](http://bio-bwa.sourceforge.net/) to map the data to
a given reference genome. The reference you use will differ depending on the species your data came from and the resources you want to use with it. For example
for a new research project consisting of Human data you would probably use the Genome Reference Consortium's [build 38 analysis set](ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines). Note that with BWA 0.7.10, mapping to alternative haplotypes has been deemed unready for production use, so you will probably wish to use
the analysis set that does not contain them.


To prepare the reference for mapping you must first index it by typing the following command where `<ref.fa>` is the path to your reference file:

	bwa index <ref.fa>

This may take several hours as it prepares the Burrows Wheeler Transform index for the reference, allowing the aligner to locate where your reads map within that reference.

Once you have finished preparing your indexed reference you can map your reads to the reference:

	bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' <ref.fa> <read1.fa> <read1.fa> > lane.sam

Typically your reads will be supplied to you in two files written in the FASTQ format.  It is particularly important to ensure that the @RG information here is correct as this information is used
by later tools. The SM field must be set to the name of the sample being processed, and LB field to the library.  The resulting mapped reads will be delivered to you in a mapping format known as [SAM](http://samtools.github.io/hts-specs/).

Because BWA can sometimes leave unusual FLAG information on SAM records, it is helpful when working with many tools to first clean up read pairing information and flags:

	samtools fixmate -O bam <lane.sam> <lane_fixmate.bam>

To sort them from name order into coordinate order:

	samtools sort -O bam -o <lane_sorted.bam> -T </tmp/lane_temp> <lane_fixmate.sam>

### Improvement
In order to reduce the number of miscalls of INDELs in your data it is helpful to realign your raw gapped alignment with the Broad's [GATK](https://www.broadinstitute.org/gatk/) Realigner.

	java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R <ref.fa> -I <lane.bam> -o <lane.intervals> --known <bundle/b38/Mills1000G.b38.vcf>
	java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R <ref.fa> -I <lane.bam> -targetIntervals <lane.intervals> --known <bundle/b38/Mills1000G.b38.vcf> -o <lane_realigned.bam>


BQSR from the Broad's GATK allows you to reduce the effects of analysis artefacts produced by your sequencing machines.  It does this in two steps, the first analyses your data to
detect covariates and the second compensates for those covariates by adjusting quality scores.

	java -Xmx4g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R <ref.fa> -knownSites >bundle/b38/dbsnp_142.b38.vcf> -I <lane.bam> -o <lane_recal.table>
	java -Xmx2g -jar GenomeAnalysisTK.jar -T PrintReads -R <ref.fa> -I <lane.bam> --BSQR <lane_recal.table> -o <lane_recal.bam>

It is helpful at this point to compile all of the reads from each library together into one BAM, which can be done at the same time as marking PCR and optical duplicates. To identify duplicates we
currently recommend the use of either the [Picard](http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates) or [biobambam's](https://github.com/gt1/biobambam) mark duplicates tool.

	java -Xmx2g -jar MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT INPUT=<lane_1.bam> INPUT=<lane_2.bam> INPUT=<lane_3.bam> OUTPUT=<library.bam>

Once this is done you can perform another merge step to produce your sample BAM files.

	samtools merge <sample.bam> <library1.bam> <library2.bam> <library3.bam>
	samtools index <sample.bam>

If you have the computational time and resources available it is helpful to realign your INDELS again:

	java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R <ref.fa> -I <sample.bam> -o <sample.intervals> --known >bundle/b38/Mills1000G.b38.vcf>
	java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R <ref.fa> -I <sample.bam> -targetIntervals <sample.intervals> --known >bundle/b38/Mills1000G.b38.vcf> -o <sample_realigned.bam>

Lastly we index our BAM using samtools:

	samtools index <sample_realigned.bam>


### Variant Calling
To convert your BAM file into genomic positions we first use mpileup to produce a BCF file that contains all of the locations in the genome.  We use this information to call genotypes
and reduce our list of sites to those found to be variant by passing this file into bcftools call.

You can do this using a pipe as shown here:

	bcftools mpileup -Ou -f <ref.fa> <sample1.bam> <sample2.bam> <sample3.bam> | bcftools call -vmO z -o <study.vcf.gz>

Alternatively if you need to see why a specific site was not called by examining the BCF, or wish to spread the load slightly you can break it down into two steps as follows:

	bcftools mpileup -Ob -o <study.bcf> -f <ref.fa> <sample1.bam> <sample2.bam> <sample3.bam>
	bcftools call -vmO z -o <study.vcf.gz> <study.bcf>

To prepare our VCF for querying we next index it using tabix:

	tabix -p vcf <study.vcf.gz>

Additionally you may find it helpful to prepare graphs and statistics to assist you in filtering your variants:

	bcftools stats -F <ref.fa> -s - <study.vcf.gz> > <study.vcf.gz.stats>
	mkdir plots
	plot-vcfstats -p plots/ <study.vcf.gz.stats>

Finally you will probably need to filter your data using commands such as:

	bcftools filter -O z -o <study_filtered..vcf.gz> -s LOWQUAL -i'%QUAL>10' <study.vcf.gz>

Variant filtration is a subject worthy of an article in itself and the exact filters you will need to use will depend on the purpose of your study and quality and depth of the data used to call the variants.

### References
* [The 1000 Genomes Project Consortium - An Integrated map of genetic variation from 1092 human genomes Nature 491, 56–65 (01 November 2012) doi:10.1038/nature11632](http://dx.doi.org/10.1038/nature11632)
* [GATK Best Practices](http://www.broadinstitute.org/gatk/guide/best-practices)

## <a name="variant_filtering"></a>Filtering of VCF files

Bcftools produces a number of parameters which may be useful for
filtering variant calls.  The most obvious of these would be the QUAL
field.

Bcftools can filter-in or filter-out using options `-i` and `-e`
respectively on the `bcftools view` or `bcftools filter1 commands.  For
example:

    bcftools filter -O z -o filtered.vcf.gz -i '%QUAL>50' in.vcf.gz
    bcftools view -O z -o filtered.vcf.gz -e 'QUAL<=50' in.vcf.gz

The quality field is the most obvious filtering method.  This is one
of the primary columns in the VCF file and is filtered using `QUAL`.
However the INFO and FORMAT fields contain many other statistics which
may be useful in distinguish true from false variants, and this is
where more complex filtering rules come in.

It can be tricky to work out the impact of various filtering rules,
and paramters may need to be changed by depth or sequencing strategy,
both technology and WGS vs Exome.  Different filtering will be needed
for SNPs and indels too.

However one useful technique, if you have a truth set available, is to
use `bcftools isec` on a VCF call file and a VCF truth file.  This
will produce 4 files containing the variants only in file 1, only in
file 2, and the variants matching in both (with the records from file 1
and in file 2 respectively).  Combining this with bcftools query will
permit construction of histograms, indicating what filtering
thresholds are appropriate.

The following are examples produced from the GIAB HG002 Illumina data
set, aligned by Novoalign.

Firstly we need to ensure both truth set and call set are normalised
using the same tool.  For this `bcftools norm -m -both -f $ref` may be
used.  Additionally you may wish to use something like `vt
decompose_blocksub` to separate out multi-allelic calls if you wish to
count each allele call separately.  If you have a bed file listing
valid regions to include or exclude, make sure to filter to those
regions too.

After this, ensure both files are bgzipped and indexed before running
isec.

```
bcftools isec -c both -p isec truth.vcf.gz call.vcf.gz
```

Now outdir/0001.vcf contains variants only found in truth.vcf.gz and
hence are false negatives.  Outdir/0002.vcf contains only variants
only in the call set, and are false positives.  Outdir/0003.vcf and
outdir/0004.vcf are the true variants.

### Depth

We may produce a histogram from outdir/0003.vcf (true) and
outdir/0001.vcf (false) to compare the distributions of the `DP`
(depth) field.  This may be either an INFO or a FORMAT field, but for
simplicitly we are restricting this guide to a single sample and using
INFO.

```
bcftools query -i 'TYPE="SNP"' -f '%DP\n ' isec/0001.vcf > DP_1
bcftools query -i 'TYPE="SNP"' -f '%DP\n ' isec/0003.vcf > DP_3
```

These files may be turned into histograms with a simple awk script or
whatever language you prefer:

```
awk '!/^\./ {a[$1]++;b++;if (min > $1) {min=$1}; if (max < $1) {max=$1}} END {print "total ",b >"/dev/stderr";for (i = min; i <= max; i++) {printf("%d\t%d\n", i, a[i])}}' DP_1 > DP_1_hist
#outputs "total  1846"
awk '!/^\./ {a[$1]++;b++;if (min > $1) {min=$1}; if (max < $1) {max=$1}} END {print "total ",b >"/dev/stderr";for (i = min; i <= max; i++) {printf("%d\t%d\n", i, a[i])}}' DP_3 > DP_3_hist
#outputs "total  257554"
```

(Note if you are also going to be selectively filtering for high
quality variants only, then you may wish to amend the "bcftools query"
command above to `-i 'TYPE="SNP" && QUAL >= 30'` to see how the
various metrics work in conjunction with quality filtering.)

Finally we can plot them in gnuplot, scaling by their total
array sums to see how the distributions of DP values differs for true
vs false SNP calls.

```
$ gnuplot
gnuplot> plot \
    "DP_1_hist" using 1:(10000*$2/1846)   with lines lw 2 title "False", \
    "DP_3_hist" using 1:(10000*$2/257554) with lines lw 2 title "True"
```

![15x normalised depth](images/15x_DP_nhist.png)
![150x normalised depth](images/150x_DP_nhist.png)

We see a sharp spike in depth for the true variants somewhere around
the expected average depth.  The false variants have a broader
distribution.  Looking at the deep data, it may appear that filtering
out variants in depth < 50 and depths > 250 would be an obvious win,
but note these plots have been normalised to have the same area so we
can clearly visualise them.  The disparity of sample sizes between
true and false variants is so large we don't really know how costly
that filtering may be.  If we don't normalise and instead plot using a
log scale then we see that even at the left (low) end there are more
true calls than false calls.

![150x log depth](images/150x_DP_lhist.png)

An alternative solution is to plot a log-odds plot, using log10(p/(1-p))
with p being the probability of this depth value being a correct call.
In this case positive values indicate how many more times likely
(logged) the event is true than false, and vice-versa for negative
values.  A small sample size correction has been applied and error
bars show the data at the right end is unreliable due to insufficient
data, but the overall picture is clear enough.

![15x log-odds depth](images/15x_DP_logodds.png)
![150x log-odds depth](images/150x_DP_logodds.png)

This now gives us a more useful way to visualise the suitable
filtering cutoffs, so we will use this style of plot from here on.

On the 15x figure it can be seen that true variants have a sharp peak,
in this case at around 14x coverage with very few true variants having
high depth.  In contrast the false variants have a long distribution
tail.  It is recommended to know the average depth and to filter at
somewhere around double that.  We observe here we would have a minimal
number of extra false negatives if we filted out `DP > 30`
(approximately 2 times the expected depth), but a considerable number
of false calls would be excluded.

At shallow data we may wish to filter at 2x average depth or slightly
higher (eg `DP > 35` on our 15x sample), while the deep data set it's
perhaps around 1.7x average depth (eg `DP > 250` on our 150x sample).


### Mapping Quality

In a heterozygous call with one allele matching the reference, the
distribution of mapping qualities for sequences matching REF versus
those matching ALT may differ due to reference bias.  This is to be
expected, however a large difference in these distributions may be
indicative of a false call.  We have a Mann-Whitney U test available
to compare these distributions.  These are normalised into a Z-score,
indicating the number of multiples of standard deviation above or
below the mean.  This is saved in the MQBZ INFO tag.

Normalised plots of these distributions can be seen here.

Note in bcftools 1.12 and earlier this is expressed as a probability
value, so filter rules will need to check against very small values,
such as `MQB < 1e-5`.

![15x log-odds MAPQ Bias](images/15x_MQBZ_logodds.png)
![150x log-odds MAPQ Bias](images/150x_MQBZ_logodds.png)

While there is a large overlap between the false and true
distributions, at both low and high depth there is a clear shifting
left for false variants.  Unfortunately the correct filtering offset
does also seem to be depth dependent. Filters of `MQBZ < -4` would be
appropriate for shallow data, and perhaps -9 for deep.

### Position

The position of a variant within the reads can matter.  We should
expect reads to be aligned fairly randomly, and thus variants to be
distributed randomly over the read.  Reference bias alignment
artifacts tend to be enriched for the ends of reads where a
substitution near the read end is usually preferable to an indel to
achieve optimal score (as alignments are pair-wise against the
reference rather than against the other sequences within in the
sample).

The RPBZ statistic is a Mann-Whitney U test represented as a Z-score
(the distance from the mean expressed in units of the standard
deviation) describing the difference in read position distributions
of REF and ALT calls.  As MQBZ this cannot be calculated for many
SNPS, but where possible it can help spot false calls due to reference
bias.

![15x log-odds Read Pos Bias](images/15x_RPBZ_logodds.png)
![150x log-odds Read Pos Bias](images/150x_RPBZ_logodds.png)

At shallow depth there isn't any discrimination power between false
and true variants.  It's more likely to be wrong at the extreme ends
of the distribution, but it's never more likely wrong than correct. At
deep data however the statistic becomes far more powerful.

The plots are largely symmetric so a filter of e.g. `RPBZ < -5 || RPBZ > +5` 
may work.

### Soft-clips

A multitude of sequence alignments having soft-clipped bases may be
indicative or a bad alignment, perhaps caused by reference bias
again.  The SCBZ is a Mann-Whitney U Z-score for the relative
distribution of length of soft-clip within the proximity of the
variant.

![15x log-odds Soft Clip Bias](images/15x_SCBZ_logodds.png)
![150x log-odds Soft Clip Bias](images/150x_SCBZ_logodds.png)

As with other Mann-Whitney U tests this statistical increases in power
as the depth increases.  It shows a sharp increase to the right end of
the distribution for false variants.

A filter of `SCBZ > 3` or `SCBZ > 4` would be appropriate for this data.


### Strand bias

This statistic is not enabled by default, but can be added with the
`-a FORMAT/SP` option of bcftools mpileup.

The plots below are normalises, and truncated in X.

![15x log-odds Strand Bias](images/15x_SP_logodds.png)
![150x log-odds Strand Bias](images/150x_SP_logodds.png)

Both true and false variants have a sharp decay, but the tail is
considerably longer for false variants.  The test is still quite
powerful for shallow data too.

As with some other metrics, the threshold is very depth specific,
needing `FORMAT/SP > 100` for the 150x data and `FORMAT/SP > 35` for
the 15x data.

### Putting it all together

While each test does not have a huge power to separate true from false
variants, some of the indicators may not be strongly correlated so
combining them together in a single clause can give a significant
boost.  If we are filtering out things that match our patterns, then
we should combine with logical OR.  For example the shallow data may use:

    bcftools view -e 'QUAL <= 10 || DP > 35 || MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || FORMAT/SP > 40 || SCBZ > 3' in.vcf

with the deep data using:

    bcftools view -e 'QUAL <= 10 || DP > 250 || MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || FORMAT/SP > 100 || SCBZ > 6' in.vcf

Note it's possible to construct some filtering rules that adjust these
thresholds according to the local depth of the data.  This is
challenging to optimise, but an example could be:

    bcftools view -e "QUAL < $qual || DP>2*$DP || MQBZ < -(3.5+4*DP/QUAL) || RPBZ > (3+3*DP/QUAL) || RPBZ < -(3+3*DP/QUAL) || FORMAT/SP > (40+DP/2) || SCBZ > (2.5+DP/30)"

Where `$qual` is the desired quality threshold and `$DP` is the
average sequencing depth.

Finally we can visualise the overall impact of our filtering by
plotting at different QUAL thresholds from 0 to 200 in increments of
10 and displaying false positives vs false negatives, with and without
the other filtering elements.  Obviously this needs a good quality
known truth set to be able to distinguish the variants.  Doing this we
see the combined power of the additional statistics.  The below figure
is for a 60x subsampling of GIAB HG002 chr1, showing SNP counts only.

![60x plot of FN vs FP](images/60x-filt.png)

The ideal position in this plot is the bottom left hand corner, with
as few false positives (high precision) and few false negatives (high
recall) as possible.  We see the green filtered plot mainly moves
leftwards, indicating a significant reduction in false positives with
a very marginal change in false negatives.


## <a name="mapping_to_cram"></a>Using CRAM within Samtools
CRAM is primarily a reference-based compressed format, meaning that only differences between the stored sequences and the reference are stored.

For a workflow this has a few fundamental effects:

1. Alignments should be kept in chromosome/position sort order. 
2. The reference must be available at all times. Losing it may be equivalent to losing all your read sequences.

Technically CRAM can work with other orders but it can become inefficient due to a large amount of random access across the reference genome. The current implementation of CRAM in htslib 1.0 is also inefficient in size for unsorted data, although this will be rectified in upcoming releases.

In CRAM format the reference sequence is linked to by the md5sum (M5 auxiliary tag) in the CRAM header (@SQ tags). This is mandatory and part of the CRAM specification. In SAM/BAM format, these M5 tags are optional. Therefore converting from SAM/BAM to CRAM requires some additional overhead to link the CRAM to the correct reference sequence.

### A Worked Example

#### Obtain some public data
We will use the first 100,000 read-pairs from a yeast data set.

	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/SRR507778/SRR507778_1.fastq.gz | gzip -d | head -100000 > y1.fastq
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/SRR507778/SRR507778_2.fastq.gz | gzip -d | head -100000 > y2.fastq
	curl ftp://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz | gzip -d > yeast.fasta

#### Prepare the BWA indices
We need to ensure there exists a .fai fasta index and also indices for whichever aligner we are using (Bwa-mem in this example).

	samtools faidx yeast.fasta
	bwa index yeast.fasta
	
#### Produce the alignments
The aligner is likely to output SAM in the same order or similar order to the input fastq files. It won't be outputting in chromosome position order, so the output is typically not well suited to CRAM.

	bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' yeast.fasta y1.fastq y2.fastq > yeast.sam

The -R option adds a read-group line and applies that read-group to all aligned sequence records. It is not necessary, but a recommended practice.

#### Sort into chromosome/positon order
Ideally at this point we would be outputting CRAM directly, but at present samtools 1.0 does not have a way to indicate the reference on the command line. We can output to BAM instead and convert (below), or modify the SAM @SQ header to include MD5 sums in the M5: field.

	samtools sort -O bam -T /tmp -l 0 -o yeast.bam yeast.sam

The "-l 0" indicates to use no compression in the BAM file, as it is transitory and will be replaced by CRAM soon. We may wish to use -l 1 if disk space is short and we wish to reduce temporary file size. 

#### Convert to CRAM format
	samtools view -T yeast.fasta -C -o yeast.cram yeast.bam

Note that since the BAM file did not have M5 tags for the reference sequences, they are computed by Samtools and added to the CRAM. In a production environment, this step can be avoided by ensuring that the M5 tags are already in the SAM/BAM header.

The last 3 steps can be combined into a pipeline to reduce disk I/O:

	bwa mem yeast.fasta y1.fastq y2.fastq | \
	samtools sort -O bam -l 0 -T /tmp - | \
	samtools view -T yeast.fasta -C -o yeast.cram -

#### Viewing in alignment and pileup format
See the variant calling workflow for more advanced examples.

	samtools view yeast.cram
	samtools mpileup -f yeast.fasta yeast.cram

### The REF_PATH and REF_CACHE

One of the key concepts in CRAM is that it is uses reference based compression. This means that Samtools needs the reference genome sequence in order to decode a CRAM file. Samtools uses the MD5 sum of the each reference sequence as the key to link a CRAM file to the reference genome used to generate it. By default Samtools checks the reference MD5 sums (@SQ "M5" auxiliary tag) in the directory pointed to by $REF_PATH environment variable (if it exists), falling back to querying the European Bioinformatics Institute (EBI) reference genome server, and further falling back to the @SQ "UR" field if these are not found.

While the EBI have an MD5 reference server for downloading reference sequences over http, we recommend use of a local MD5 cache. We have provided with Samtools a basic script (misc/seq_cache_populate.pl) to convert your local yeast.fasta to a directory tree of reference sequence MD5 sums:

	<samtools_src_dir>/misc/seq_cache_populate.pl -root /some_dir/cache yeast.fasta
	export REF_PATH=/some_dir/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
	export REF_CACHE=/some_dir/cache/%2s/%2s/%s

REF_PATH is a colon separated list of directories in which to search for files named after the sequence M5 field. The : in http:// is not considered to be a separator. Hence using the above setting, any CRAM files that are not cached locally may still be looked up remotely.

In this example "%2s/%2s/%s" means the first two digits of the M5 field followed by slash, the next two digits and slash, and then the remaining 28 digits.  This helps to avoid one large directory with thousands of files in it.

The REF_CACHE environment variable is used to indicate that any downloaded reference sequences should be stored locally in this directory in order to avoid subsequent downloads.  This should normally be set to the same location as the first directory in REF_PATH.
