---
layout: default
title: Samtools - Documentation
---
## Manual pages

Documentation for BCFtools, SAMtools, and HTSlib's utilities is available
by using <code>man <em>command</em></code> on the command line.
The manual pages for several releases are also included below --- be sure
to consult the documentation for the release you are using.

* [bcftools 1.3.1](bcftools.html) (older versions:
      [1.3](bcftools-1.3.html),
      [1.2](bcftools-1.2.html),
      [1.1](bcftools-1.1.html),
      [1.0](bcftools-1.0.html),
      [0.1.19](samtools-0.1.19.html "included in samtools-0.1.19"))
* [bgzip 1.3.1](tabix.html) (older versions:
      [1.3](tabix-1.3.html),
      [1.2.1](tabix-1.2.1.html),
      [1.1](tabix-1.1.html))
* [htsfile 1.3.1](htsfile.html) (older version:
      [1.3](htsfile-1.3.html),
      [1.2.1](htsfile-1.2.1.html))
* [samtools 1.3.1](samtools.html) (older versions:
      [1.3](samtools-1.3.html),
      [1.2](samtools-1.2.html),
      [1.1](samtools-1.1.html),
      [1.0](samtools-1.0.html),
      [0.1.19](samtools-0.1.19.html))
* [tabix 1.3.1](tabix.html) (older versions:
      [1.3](tabix-1.3.html),
      [1.2.1](tabix-1.2.1.html),
      [1.1](tabix-1.1.html),
      [1.0](tabix-1.0.html))

## File formats

SAMtools conforms to the specifications produced by the GA4GH File Formats working group. Details of the current specifications are available on the  [hts-specs page](http://samtools.github.io/hts-specs).

HTSlib also includes brief manual pages outlining aspects of several of
the more important file formats.
These are available via <code>man <em>format</em></code> on the command line
or here on the web site:

* [faidx](faidx.html) describes _.fai_ FASTA index files
* [sam](sam.html) lists the mandatory SAM fields and meanings of flag values
* [vcf](vcf.html) lists the mandatory VCF fields and common INFO tags

## Benchmarks

* [Zlib implementations](../benchmarks/zlib.html) comparing samtools read and
  write speeds.

* [CRAM comparisons](../benchmarks/CRAM.html) between version 2.1,
  version 3.0 and BAM formats.

## Publications

* Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. _[PMID: [19505943](http://www.ncbi.nlm.nih.gov/pubmed/19505943)]_
* Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. _[PMID: [21903627](http://www.ncbi.nlm.nih.gov/pubmed/21903627)]_
* Danecek P., Schiffels S., Durbin R. Multiallelic calling model in bcftools (-m) _[[link](http://samtools.github.io/bcftools/call-m.pdf)]_
* Li H. Improving SNP discovery by base alignment quality. Bioinformatics. 2011 Apr 15;27(8):1157-8. doi: 10.1093/bioinformatics/btr076. Epub 2011 Feb 13. _[PMID: [21320865](http://www.ncbi.nlm.nih.gov/pubmed/21320865)]_
* Durbin R. Segregation based metric for variant call QC _[[link](http://samtools.github.io/bcftools/rd-SegBias.pdf)]_
* Li H, Mathematical Notes on SAMtools Algorithms _[[link](http://www.broadinstitute.org/gatk/media/docs/Samtools.pdf)]_
