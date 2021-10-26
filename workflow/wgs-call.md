---
permalink: /workflow/wgs-call.html
layout: default
---

# WGS/WES Mapping to Variant Calls

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

