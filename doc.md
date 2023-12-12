---
layout: default
title: Samtools - Documentation
---
## Manual pages

Documentation for BCFtools, SAMtools, and HTSlib's utilities is available
by using <code>man <em>command</em></code> on the command line.
The manual pages for several releases are also included below --- be sure
to consult the documentation for the release you are using.

* [annot-tsv](annot-tsv.html)
* [bcftools](bcftools.html)
* [bgzip](bgzip.html)
* [htsfile](htsfile.html)
* [samtools](samtools.html)
* [tabix](tabix.html)

Older manual pages are available for releases: [0.1.19](0.1.19), [1.0](1.0), [1.1](1.1), [1.2](1.2), [1.3](1.3), [1.3.1](1.3.1), [1.4](1.4), [1.4.1](1.4.1), [1.5](1.5), [1.6](1.6), [1.7](1.7), [1.8](1.8), [1.9](1.9), [1.10](1.10), [1.11](1.11), [1.12](1.12), [1.13](1.13), [1.14](1.14), [1.15](1.15), [1.16](1.16), [1.17](1.17), [1.18](1.18)

## HowTos

* Options for [viewing files with & without headers](../howtos/headers.html)

### [BCFtools](http://samtools.github.io/bcftools/)

- Installation
    - [Latest development version](http://samtools.github.io/bcftools/howtos/install.html)
- Calling
    - [CNV calling](http://samtools.github.io/bcftools/howtos/cnv-calling.html)
    - [Consequence calling](http://samtools.github.io/bcftools/howtos/csq-calling.html)
    - [Consensus calling](http://samtools.github.io/bcftools/howtos/consensus-sequence.html)
    - [ROH calling](http://samtools.github.io/bcftools/howtos/roh-calling.html)
    - [Variant calling and filtering](http://samtools.github.io/bcftools/howtos/variant-calling.html)

- Tips and Tricks
    - [Converting formats](http://samtools.github.io/bcftools/howtos/convert.html)
    - [Extracting information](http://samtools.github.io/bcftools/howtos/query.html)
    - [Filtering expressions](http://samtools.github.io/bcftools/howtos/filtering.html)
    - [Plugins](http://samtools.github.io/bcftools/howtos/plugins.html)




## File formats

SAMtools conforms to the specifications produced by the GA4GH File Formats working group. Details of the current specifications are available on the  [hts-specs page](http://samtools.github.io/hts-specs).

HTSlib also includes brief manual pages outlining aspects of several of
the more important file formats.
These are available via <code>man <em>format</em></code> on the command line
or here on the web site:

* [faidx](faidx.html) describes _.fai_ FASTA index files
* [sam](sam.html) lists the mandatory SAM fields and meanings of flag values
* [vcf](vcf.html) lists the mandatory VCF fields and common INFO tags
* [htslib-s3-plugin](htslib-s3-plugin.html) describes the S3 plugin

## Algorithms

* [Duplicate](../algorithms/duplicate.html) marking.

## Benchmarks

* [Zlib implementations](../benchmarks/zlib.html) comparing samtools read and
  write speeds.

* [CRAM comparisons](../benchmarks/CRAM.html) between version 2.1,
  version 3.0 and BAM formats.

## Publications

### Software Packages

A joint publication of SAMtools and BCFtools improvements over the last 12
years was published in 2021.

* Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H, **Twelve years of SAMtools and BCFtools**, *GigaScience* (2021) 10(2) giab008 [[33590861](https://pubmed.ncbi.nlm.nih.gov/33590861)]

The same journal issue also saw the HTSlib paper, describing the C library.

* Bonfield JK, Marshall J, Danecek P, Li H, Ohan V, Whitwham A, Keane T, Davies RM, **HTSlib: C library for reading/writing high-throughput sequencing data**, *GigaScience* (2021) 10(2) giab007 [[33594436](https://pubmed.ncbi.nlm.nih.gov/33594436)]


### File formats

The introduction of the **SAM/BAM format** and the `samtools` command line tool:

* Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, **The Sequence alignment/map (SAM) format and SAMtools**, *Bioinformatics* (2009) 25(16) 2078-9 [[19505943](http://www.ncbi.nlm.nih.gov/pubmed/19505943)]

Extension of the SAM/BAM format to support *de novo* assemblies:

* Cock PJA, Bonfield JK, Chevreux B, Li H, **SAM/BAM format v1.5 extensions for *de novo* assemblies**, *bioRxiv* (2015) 020024 [[doi:10.1101/020024](http://dx.doi.org/10.1101/020024)]

The introduction of the **CRAM format**:

* Hsi-Yang Fritz M, Leinonen R, Cochrane G, and Birney E, **Efficient storage of high throughput DNA sequencing data using reference-based compression**, *Genome Research* (2011) 21(5) 734-740. [[21245279](http://www.ncbi.nlm.nih.gov/pubmed/21245279)]

The introduction of the **VCF format**:

* Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R, 1000 Genomes Project Analysis Group, **The variant call format and VCFtools**, *Bioinformatics* (2011) 27(15) 2156-8 [[21653522](http://www.ncbi.nlm.nih.gov/pubmed/21653522)]

### Calling and analysis

The original mpileup calling algorithm plus mathematical notes (`mpileup/bcftools call -c`):

* Li H, **A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data**, *Bioinformatics* (2011) 27(21) 2987-93. [[21903627](http://www.ncbi.nlm.nih.gov/pubmed/21903627)]
* Li H, **Mathematical Notes on SAMtools Algorithms** (2010) [[link](http://www.broadinstitute.org/gatk/media/docs/Samtools.pdf)]

Mathematical notes for the updated multiallelic calling model (`mpileup/bcftools call -m`):

* Danecek P, Schiffels S, and Durbin R, **Multiallelic calling model in bcftools (-m)** (2014) [[link](http://samtools.github.io/bcftools/call-m.pdf)]

Hidden Markov model for detecting runs of homozygosity (`bcftools roh`):

* Narasimhan V, Danecek P, Scally A, Xue Y, Tyler-Smith C, and Durbin R, **BCFtools/RoH: a hidden Markov model approach for detecting autozygosity from next-generation sequencing data**, *Bioinformatics* (2016) 32(11) 1749-51 [[26826718](http://www.ncbi.nlm.nih.gov/pubmed/26826718)]

Copy number variation/aneuploidy calling from microarray data (`bcftools cnv/bcftools polysomy`):

* Danecek P, McCarthy SA, HipSci Consortium, and Durbin R, **A Method for Checking Genomic Integrity in Cultured Cell Lines from SNP Genotyping Data**, *PLoS One* (2016) 11(5) e0155014 [[27176002](http://www.ncbi.nlm.nih.gov/pubmed/27176002)]

Haplotype-aware calling of variant consequences (`bcftools csq`):

* Danecek P, McCarthy SA, **BCFtools/csq: Haplotype-aware variant consequences**, *Bioinformatics* (2017) 33(13) 2037-39 [[28205675](http://www.ncbi.nlm.nih.gov/pubmed/28205675)]

### Other

Base alignment quality (BAQ) method improve SNP calling around INDELs:

* Li H, **Improving SNP discovery by base alignment quality**, *Bioinformatics* (2011) 27(8) 1157-8 [[21320865](http://www.ncbi.nlm.nih.gov/pubmed/21320865)]

Segregation based QC metric originally implemented in [SGA](https://github.com/jts/sga):

* Durbin R, **Segregation based metric for variant call QC** (2014) [[link](http://samtools.github.io/bcftools/rd-SegBias.pdf)]
