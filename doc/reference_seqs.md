---
permalink: /doc/reference_seqs.html
layout: default
title: HTSlib - Reference Sequences
highlighting: yes
---

# Reference Sequences

## Introduction

CRAM typically uses reference based compression where individual bases in aligned records are compared against a known reference sequence and only bases differing are stored.
This is used as a form of sequence compression, but it can also help construct the `MD` and `NM` auxiliary tags, further saving storing costs.

Most CRAMs will use an external reference, however it is also possible to use an internal embedded reference within the CRAM file or to use no reference at all.

## Search Order For Finding References

When looking for an external reference the search order used by HTSlib is as follows.  

1. Use any reference FASTA file specified by the `hts_set_fai_filename()` function or samtools command line options, such as `samtools view -T`.

2. Search the `REF_CACHE` directory for the sequence MD5sum listed in the `@SQ` header.  This data will be unformatted raw sequence.

3. In turn, search each `REF_PATH` location element for the sequence MD5sum.  This can be either local directories or remote URLs.

4. Look for a local FASTA file listed in the `@SQ` header UR tag.  Earlier versions of HTSlib unintentionally also utilised http and ftp based UR tags.  These are no longer supported because they do not operate in conjunction with REF_CACHE.

5. Finally if none of the above succeed, the library will return a failure code.

The sections below explain these variables and formats in more detail.


## External References

An external reference is one where the reference is not stored as part of the CRAM file.
This offers optimal compression, provided the sample genome matches the external reference sufficiently well.
It is required to be known for both compression and decompression.

When using an external reference it is important for the software to be able to find the appropriate reference file.
We can specify this manually using e.g. `samtools view -T reference.fa`.
Specifying the incorrect reference would lead to incorrect sequence decompression, but each CRAM slice also includes the MD5sum of the portion of reference used during encode and this is validated during decode.
Hence supplying the wrong reference will typically lead to an "MD5sum checksum reference mismatch" error.

A better solution is to find references by their whole-chromosome MD5sum, as recorded in the `@SQ M5:` SAM header field.
The `REF_PATH` environment variable may be used as a colon-separated (semi-colon on Windows) list of locations to search for this MD5sum.
The data returned is expected to be the unadulterated upper-cased sequence with no extra formatting or line-wrapping.
The locations in `REF_PATH` may be directories or a URL to a GA4GH refget server.
If using a filesystem directory path, to avoid thousands or even millions of reference sequences all in the same directory it is advisable to separate the MD5sum filename into a series of subdirectories.
This is achieved by using %*num*s within REF_PATH where *num* is the next number of digits from the MD5sum string.
For example to find GRCh38 chr1 with an MD5sum of 6aef897c3d6ff0c78aff06ac189178dd and having `REF_PATH=/path/to/references/%2s/%2s/%s`, HTSlib will look for filename `/path/to/references/6a/ef/897c3d6ff0c78aff06ac18917`.

We can use URLs as part of the `REF_PATH` search path.
Given these naturally have colons in `https://` and perhaps for a port number, these are considered to be a single element when tokenising the path.
However this is limited to http, https and ftp URIs so it is recommended to escape literal colons within URIs by using a double-colon instead.

Note it is faster to fetch sequences from an MD5sum-structured local directory than fetching it from local FASTA file as there is no parsing or validation required on the sequences.
Indeed HTSlib typically accesses the MD5sum sequence using the UNIX `mmap()` system call, which also reduces memory usage when many processes are accessing the same reference file.

## Reference Caching

Fetching references from a remote refget server can be slow and potentially puts a significant burden on the external resource.
It is recommended that we use a combination of a local directory-based reference directory structure prior and a locally-caching refget server.
For example `REF_PATH=/local/refs/%2s/%2s/%2s:https://refget.server/sequence::8000/%s`.

HTSlib comes with its own ref-cache tool to act as a caching proxy to the EBI's CRAM registry (or any other refget server specified).
See the ref-cache.1 man page for further details.

## Prepopulating a reference cache

If the `REF_CACHE` environment variable is set, any non-local URIs used for fetching a reference via `REF_PATH` will be saved to the location specified in `REF_CACHE`.
This environment variable shares the same percent-s encoding rules as `REF_PATH`.
Note `REF_CACHE` is not appropriate for a resource shared between multiple users as this may lead to permission problems.

We can pre-populate our /local/refs directory from a set of FASTA files using Samtools' seq_cache_populate.pl script.
This can be particularly useful when processing data in cloud environments where we wish to have no external network access beyond the initial provisioning of the cloud environment.
For example:

```
seq_cache_populate.pl -root=/local/refs -subdirs 2 *.fa
```

At the Wellcome Sanger Institute we use a local web server as the final REF_PATH element.
This web server proxies the EBI refget server while also adding the entry to the local directory-based REF_PATH location, thus ensuring our local cache is updated and avoiding subsequent external access.

## Internal / Embedded References

CRAM can store the reference internally as part of the file.
This means the reference will not need to be known during decompression.
While this feature was added for efficient storage of denovo assemblies, it may still be worth considering for data aligned against widely known reference genomes when wanting a self contained CRAM that has no requirements on external data.

This may seem to remove the benefits of reference based encoding, but there will usually be many alignments covering the same location so this is a form of deduplication analogous to LZ compression techniques.
Each CRAM slice stores the portion of the reference covered by alignments within that slice.
This obviously makes internal / embedded references only compatible with coordinate sorted aligned records.
The reference stored in the slice can be directly copied from the reference used during alignment, or created on-the-fly by computing the consensus from the alignments.
This is controlled in HTSlib with the `embed_ref` or `embed_ref=1` option, which uses an externally supplied reference, and `embed_ref=2` which creates a new consensus reference.
The latter typically gives better compression of sequence data but may harm MD and NM compression as these describe the delta against the reference used during alignment.
However on-the-fly consensus reference has the benefit of working even when we do not know which reference was used during BAM alignment.
That is not an ideal situation to be in, but if we cannot find the external reference then the CRAM encoder will automatically switch to building its own consensus instead.

The space taken up by the reference in a CRAM file is typically around 2 bits per covered base, irrespective of alignment depth.
This means a CRAM file for a human genome using an embedded reference will be around 700MB larger than the same CRAM using an external reference.

For example `samtools view -O cram,embed_ref=2 in.bam -o out.cram` will create a CRAM file using a embedded consensus as the reference.


## Benchmarks

Cram sizes for different reference modes are listed below.
This data file is 10 million NovaSeq X records aligned the first 47Mbp of human chromosome 1, at approximately 50x coverage.

{:.table}
Size (MB) | Format         | Options
--------- | -------------- | -------
515.1     | BAM            |
175.5     | CRAM 3.1       | (default)
186.8     | CRAM 3.1       | embed_ref=1
185.8     | CRAM 3.1       | embed_ref=2
224.2     | CRAM 3.1       | no_ref

