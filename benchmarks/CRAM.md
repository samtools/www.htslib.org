---
permalink: /benchmarks/CRAM.html
layout: default
title: CRAM Benchmarks
highlighting: yes
---
CRAM benchmarking
=================

(Updated Sep 2019)

<!-- Tested on my gen1c (Intel(R) Xeon(R) CPU E5-2660 0 @ 2.20GHz)

R/W:  perf stat ./test/test_view -@8 -C -o opt... in.bam -p out.cram
READ: perf stat ./test/test_view -@8 -B out.cram
-->

Benchmarks of CRAM 2.1 and 3.0 are using the faster CRAM codecs;
primarily deflate and rANS.

Also included is the performance of the proposed CRAM v3.1 standard.
This is not yet a ratified GA4GH standard, but these figures give
indicative results.

The options listed below also include the new proposed compression
profiles (fast, normal (default), small and archive) that ease the
trade-off between speed vs size vs random access.  The profiles are
synonyms for a collection of existing options.  At the time of
writing, these profiles are:

{:.table}
Profile          | CRAM versions | options
---------------- | ------------- | -------
fast             | 3.0, 3.1      | seqs_per_slice=1000,  level=1
normal (default) | 3.0, 3.1      | seqs_per_slice=10000
small            | 3.0           | seqs_per_slice=25000, level=6,use_bzip2
small            | 3.1           | seqs_per_slice=25000, level=6,use_bzip2,use_fqz
archive          | 3.0           | seqs_per_slice=100000,level=7,use_bzip2
archive          | 3.1           | seqs_per_slice=100000,level=7,use_bzip2,use_fqz,use_arith

To demonstrate the absolute smallest size we use add option "use_lzma"
to the archive profile tests.  This adds considerable encode cost, but
minimal decode.


Coordinate sorted human data
----------------------------

This test set is chr1 of NA12878_S1, downloaded from
<http://www.ebi.ac.uk/ena/data/view/ERP002490>

Conversion from BAM to <format> is achieved via htslib test_view -b or
-C using 8 threads on a multi-core Intel Xeon E5-2660 system.
Decoding times are computed using test_view -B (benchmarking mode)
with performs input, uncompression and decoding only.  Encoding
time includes decoding of the input BAM, so subtract the BAM decoding
time to get encode-only, although the difference is only minor.

All times are reported as wall-clock, although some I/O time will
impact this test as the file is large, particularly decode times.
In this first test we also show both vanilla Zlib and Libdeflate
implementations of the gzip standard.


{:.table}
Format               | Options          |     Size(Mb) | Encoding(s) | Decoding(s)
-------------------- | ---------------- | -----------: | ----------: | ----------:
BAM (zlib)           |                  |       121710 |        3101 |         438
BAM (libdeflate)     |                  |       122405 |        2192 |         357
-------------------- | ---------------- | -----------: | ----------: | ----------:
CRAM v2.1            |                  |        78259 |        2121 |         407
-------------------- | ---------------- | -----------: | ----------: | ----------:
CRAM v3.0            |                  |        66494 |        1094 |         407
CRAM v3.0            | small            |        64682 |        2122 |         573
CRAM v3.0            | archive,use_lzma |        63629 |        6654 |         394
-------------------- | ---------------- | -----------: | ----------: | ----------:
CRAM v3.1 (proposed) |                  |        62150 |        1324 |         417
CRAM v3.1 (proposed) | small            |        56204 |        2465 |        1405
CRAM v3.1 (proposed) | archive,use_lzma |        54237 |        3048 |        1395

E.Coli: sorted, unsorted, with and without references
-----------------------------------------------------

To compare CRAM efficiency in a variety of circumstances we chose a
smaller dataset to more completely explore the parameter space.
MiSeq_Ecoli_DH10B_110721_PF.bam is the smallest example data taken
from the Deez paper, so we also include Deez itself here for
comparison too.

The BAM/SAM.gz implementation here uses libdeflate.  Again we use 8
threads, but note deez was only able to use around 2.

### Position sorted

{:.table}
Format       | Options          |   Size(Mb) | Encoding(s) | Decoding(s)
------------ | ---------------- | ---------: | ----------: | ----------:
BAM          |                  |       1420 |        20.4 |        1.8
SAM.gz       |                  |       1387 |        22.4 |        2.9
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v2.1    |                  |       1048 |        24.6 |        3.5
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v3.0    |                  |        868 |        10.6 |        3.8
CRAM v3.0    | small            |        844 |        23.6 |        4.8
CRAM v3.0    | archive,use_lzma |        838 |        94.9 |        4.5
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v3.1    |                  |        830 |        13.2 |        4.2
CRAM v3.1    | small            |        771 |        35.4 |       17.3
CRAM v3.1    | archive,use_lzma |        747 |        81.7 |       16.9
------------ | ---------------- | ---------: | ----------: | ----------:
Deez         |                  |        842 |       134.6 |       62.8

With light-weight level 1 compression and uncompressed level 0 files
we see CRAM 3 being slower for uncompressed data than CRAM 2.  This is
due to the additional CRC integrity checks.

Note the slower time for BAM level 0 than level 1 is purely down to
increased disk I/O costs; CPU times double for level 1.  Why SAM does
not pay this penalty is unknown, but it is likely this picture would
change given a large enough file.

{:.table}
Format    | Options |   Size(Mb) | Encoding(s) | Decoding(s)
--------- | ------- | ---------: | ----------: | ----------:
SAM.gz    | level=1 |       1830 |       22.7  |        3.2
BAM       | level=1 |       1531 |       17.8  |        1.9
CRAM v2.1 | level=1 |       1126 |       12.2  |        3.5
CRAM v3.0 | level=1 |        889 |        9.5  |        3.6
CRAM v3.1 | level=1 |        858 |       10.3  |        3.7
--------- | ------- | ---------: | ----------: | ----------:
SAM       | level=0 |       4463 |        7.8  |        1.5
BAM       | level=0 |       4310 |       19.6  |        2.1
CRAM v2.1 | level=0 |       2632 |       10.1  |        4.1
CRAM v3.0 | level=0 |       2632 |       11.4  |        4.9
CRAM v3.1 | level=0 |       2632 |       10.9  |        4.8

#### Embeded reference & referenceless encoding (position sorted alignments)

By default aligned CRAM uses an external reference file.  Portions of
that reference can be embedded within each slice to remove this
external file dependency.  On deep data this has minimal impact as the
reference is small in comparison to the alignments.

CRAM can also do non reference-based compression, storing the sequence
as-is (like BAM).  This leads to larger files.

{:.table}
Format    | Options   |   Size(Mb) | Encoding(s) | Decoding(s)
--------- | --------- | ---------: | ----------: | ----------:
CRAM v2.1 |           |       1048 |        24.6 |        3.5
CRAM v2.1 | embed_ref |       1049 |        25.2 |        3.8
CRAM v2.1 | no_ref    |       1089 |        53.7 |       13.9
--------- | --------- | ---------: | ----------: | ----------:
CRAM v3.0 |           |        868 |        10.6 |        3.8
CRAM v3.0 | embed_ref |        869 |        10.6 |        3.7
CRAM v3.0 | no_ref    |        905 |        13.4 |        3.8
--------- | --------- | ---------: | ----------: | ----------:
CRAM v3.1 |           |        830 |        13.2 |        4.2
CRAM v3.1 | embed_ref |        831 |        13.4 |        4.0
CRAM v3.1 | no_ref    |        866 |        15.9 |        4.6

The significant speed difference of no_ref between version 2.1 and 3.0
is due to improved ways of storing multi-base differences instead of
requiring one CRAM feature for each base call.

### Name sorted alignments

This is the same file above, with aligned sequencing data, but sorted
into name order using "samtools sort -n".  BAM is significantly larger
as the sequences are no longer in sorted order, harming gzip, but CRAM
does not change size considerably.  This is due to the use of
reference based compression.  With referenceless compression CRAM will
grow in size, similar to BAM, as is visible with the "no_ref" option.
Although "no_ref" it makes minimal difference here, with a very large
reference it may be preferable to use this on name sorted data to
reduce memory usage.

{:.table}
Format       | Options          |   Size(Mb) | Encoding(s) | Decoding(s)
------------ | ---------------- | ---------: | ----------: | ----------:
BAM          |                  |       2016 |        26.2 |        2.3
SAM.gz       |                  |       2005 |        30.8 |        3.4
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v2.1    |                  |       1034 |        26.3 |        5.1
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v3.0    | no_ref           |       1327 |        12.7 |        4.0
CRAM v3.0    |                  |        846 |        12.5 |        4.7
CRAM v3.0    | small            |        832 |        22.9 |        4.4
CRAM v3.0    | archive,use_lzma |        818 |        84.7 |        4.8
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v3.1    | no_ref           |       1302 |        15.5 |        4.2
CRAM v3.1    |                  |        825 |        14.6 |        5.9
CRAM v3.1    | small            |        768 |        37.4 |       18.1
CRAM v3.1    | archive,use_lzma |        741 |        98.3 |       17.5
 
 
### Unmapped (name order)

This is the name sorted data above, but with alignments and all
auxiliary tags stripped out.  This was achieved by converting back to
FASTQ via "samtools fastq" and from there back to unaligned BAM.  As
expected the CRAMs are broadly similar in size to the no_ref mapped
name sorted files.

{:.table}
Format       | Options          |   Size(Mb) | Encoding(s) | Decoding(s)
------------ | ---------------- | ---------: | ----------: | ----------:
FASTQ.pigz   |                  |       1681 |        58.0 |       21.4
FASTQ.bgzf   |                  |       1710 |        30.0 |        7.9
BAM          |                  |       1686 |        21.7 |        2.0
SAM.gz       |                  |       1722 |        28.3 |        2.6
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v2.1    |                  |       1476 |        28.3 |        3.6
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v3.0    |                  |       1236 |        12.0 |        3.7
CRAM v3.0    | small            |       1231 |        19.5 |        3.5
CRAM v3.0    | archive,use_lzma |       1041 |       378.5 |        5.2
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v3.1    |                  |       1208 |        14.6 |        5.9
CRAM v3.1    | small            |       1139 |        73.0 |       16.8
CRAM v3.1    | archive,use_lzma |        967 |       423.3 |       18.2
 
Note the fastq was compressed with pigz and bgzip both using 8
threads.  Pigz is smaller, but bgzip's use of libdeflate greatly
improves the speed.  Both are significantly behind unaligned CRAM
though so our recommendation is against storing data in native FASTQ
format.

### Unmapped (Minhash order)

As above, but passed through the experimental "samtools sort -M"
command first.  This clusters reads by a hash of their sequence,
having the effect of grouping similar looking data together which
helps LZ compression algorithms.  Note the data is still unaligned.
It is a quick alternative to (a better) full genome assembly.  With 8
threads this sort process took 21 seconds real time, 170 seconds CPU,
although expect this to be less performant on a very large file as
would spill temporary files to disk and require a large merge sort.

Note some aligners will need these files sorting (or collating) back
to name order prior to converting back to FASTQ.

{:.table}
Format       | Options          |   Size(Mb) | Encoding(s) | Decoding(s)
------------ | ---------------- | ---------: | ----------: | ----------:
BAM          |                  |       1207 |        17.0 |        1.5
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v3.0    |                  |        871 |        13.3 |        3.6
CRAM v3.0    | archive,use_lzma |        831 |       125.8 |        4.9
------------ | ---------------- | ---------: | ----------: | ----------:
CRAM v3.1    |                  |        842 |        16.3 |        3.6
CRAM v3.1    | archive,use_lzma |        748 |       163.5 |       16.0

The effect of "sort -M" on a small deeply sequenced genome is
profound, giving file sizes comparable to the aligned position sorted
data and around half the size of the name sorted compressed FASTQ
file.  Expect this effect to be less pronounced on shallow data sets
or much larger genomes.
