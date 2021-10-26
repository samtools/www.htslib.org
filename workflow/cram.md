---
permalink: /workflow/cram.html
layout: default
---

# Using CRAM within Samtools

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
