---
layout: default
title: Samtools - Workflows: FASTQ to BAM / CRAM
redirect_from: /workflow#fastq_to_bam
---

# <a name="fastq_to_bam"></a>FASTQ to BAM / CRAM - Version 1.0

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


