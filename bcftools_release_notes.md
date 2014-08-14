---
layout: default
title: Samtools - BCFTools 1.0 Release Notes
---

BCFtools 1.0 RELEASE NOTES
==========================

This is a major new release of bcftools, a set of utilities for variant calling 
and manipulating VCF and BCF files. The bcftools package has been separated from 
the samtools package and has undergone a major rewrite to be based on the open 
source C library [htslib].

* **BCF2**. Version 1.0 is designed to work with VCF and BCF2 (uncompressed and 
bgzip compressed). The BCF1 format output by versions of samtools <=0.1.19 is 
*not* compatible with bcftools 1.0. To read BCF1 files one can use the view 
command from old versions of bcftools packaged with samtools <=0.1.19 to convert 
to VCF, which can then be read by bcftools 1.0. Currently bcftools outputs and 
supports BCF2.2. See the [specifications] page.

* **VARIANT CALLING**: Variant calling from the output of the `mpileup` command 
in samtools (1.0) is now done using `bcftools call` rather than `bcftools view`. 
Users are now required to choose between the old samtools calling model 
(`-c/--consensus-caller`) and the new multi-allelic calling model 
(`-m/--multiallelic-caller`). The multi-allelic calling model is recommended for 
most tasks and has been improved to correctly handle priors. There is also new 
support to output gVCF.

* `bcftools view` is no longer used for variant calling. Instead, it is now used 
for VCF<=>BCF conversion, viewing, subsetting and filtering VCF and BCF files.

* The `bcftools index` command creates CSI (coordinate-sorted index) by default. 
The CSI format supports indexing of chromosomes up to length 2^31 (see the 
[specifications]). TBI (tabix index) index files, which support chromosome 
lengths up to 2^29, can still be created by using the `-t/--tbi` option or using 
the `tabix` program packaged with [htslib]. When loading an index file, bcftools 
will try the CSI first and then the TBI.

* The bcftools commands `cat`, `ld` and `ldpair` are no longer supported.

* Most of the Perl-based tools from [vcftools] have been converted to use htslib 
and are now part of bcftools and are now much faster. See the man page for 
details on using these tools.

[htslib]: http://htslib.org/

[specifications]: http://samtools.github.io/hts-specs/

[vcftools]: http://vcftools.sourceforge.net/
