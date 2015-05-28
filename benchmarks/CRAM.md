---
permalink: /benchmarks/CRAM.html
layout: default
title: CRAM Benchmarks
highlighting: yes
---
CRAM benchmarking
=================

Benchmarks are using the faster CRAM codecs; primarily deflate and
rANS.  For comparison we also include "Io_lib"'s Scramble tool for
bzip2 and lzma CRAM (not yet supported in htslib) and the Deez tool on
one data set.

Coordinate sorted human data
----------------------------

This test set is chr1 of NA12878_S1, downloaded from
<http://www.ebi.ac.uk/ena/data/view/ERP002490>

Conversion from BAM to <format> is achieved via htslib test_view -b or
-C.  Decoding times are computed using test_view -B (benchmarking
mode) with performs input, uncompression and decoding only.  Encoding
time includes decoding of the input BAM, so subtract the BAM decoding
time to get encode-only, although the difference is only minor.

All times are reported as wall-clock, although typically these
algorithms are CPU bound so the cpu time is comparable.

{:.table}
Format   |       Size | Encoding(s) | Decoding(s)
-------- | ---------: | ----------: | ----------:
BAM      | 8164228924 |        1216 |         111
CRAM v2  | 5716980996 |        1247 |         182
CRAM v3  | 4922879082 |         574 |         190

E.Coli: sorted, unsorted, with and without references
-----------------------------------------------------

To compare CRAM efficiency in a variety of circumstances we chose a
smaller dataset to more completely explore the parameter space.
MiSeq_Ecoli_DH10B_110721_PF.bam is the smallest example data taken
from the Deez paper, so we also include Deez here for comparison.

###Position sorted

{:.table}
Format       |       Size | Encoding(s) | Decoding(s) | Notes
------------ | ---------: | ----------: | ----------: | :--------------
SAM          | 5579036306 |          46 |         -   |
BAM          | 1412001095 |         209	|	17.7  |
CRAM v2      | 1053744556 |         183 |       26.9  |
CRAM v3      |  869500447 |          75 |       30.8  |
CRAM v3+bz2  |  850165878 |         124 |       45.4  | Via Scramble -j
Deez         |  870040062 |         208 |      165.0  | Via deez

Extra decoding time for CRAM v3 is largely explained by the additional
CRC checksums.

The effect of varying compression levels:

{:.table}
Format  | Level     |        Size | Encoding(s)
------- | --------- | ----------: | ----------:
BAM     | 9         |  1399448787 |         403
BAM     | (default) |  1412001095 |         209
BAM     | 1         |  1616365585 |          88
BAM     | u         |  4310414959 |          45
CRAMv3  | 9         |   862152172 |         193
CRAMv3  | (default) |   869500447 |          75
CRAMv3  | 3         |   881163788 |          68
CRAMv3  | 1         |   886716515 |          69
CRAMv3  | u         |  2631665078 |          45

Compression level "u" is uncompressed.  Note there is almost no
difference in speed between CRAM level 1 and the default level.  Maybe
we need to make -1 more aggressively fast at the expense of ratio.

Also note that BAM -1 is slower to encode than CRAMv3 at default
levels, although it will be faster to decode.  That makes me wonder
about how we should deal with temporary files. (Ideally with neither
BAM nor CRAM compression, but LZ4 or Snappy.)

###Name sorted

{:.table}
Format  | Size       | Encoding(s)
------- | ---------: | ----------:
BAM     | 2005274489 | 284
CRAMv2  | 1046736335 | 206
CRAMv3  |  848319136 |  99


###Embedding & Reference-less encoding

Embedded reference - no external file dependencies:

{:.table}
Format  | Size       |  Encoding(s)
------- | ---------: | ----------:
CRAM v2 | 1055144162 | 187
CRAM v3 |  870779925 |  78

Note that this is almost the same as the default mode of using an
external reference as the sequence depth is high.


Non-reference encoding - all sequence bases are stored verbatim:

{:.table}
Format  | Size       | Encoding(s)
------- | ---------: | ----------:
CRAM v2 | 1140863900 | 387
CRAM v3 |  941857968 | 100

The significant speed difference between version 2.1 and 3.0 is due to
improved ways of storing multi-base differences instead of requiring
one CRAM feature for each base call.


Unmapped data
-------------

Human gut sample SAMEA728920 from http://www.ebi.ac.uk/ena/data/view/ERA000116
This is unmapped data, converted from FASTQ to SAM via biobambam.


{:.table}
Format       | Size        | Encoding(s)  | Decoding(s) | Notes
------------ | ---------:  | ----------:  | ----------: | --------------------------------
SAM          | 1443985040  |  15	  |     -       | I/O bound
BAM          |  428540917  |  64	  |     5.5     | 
CRAM v2      |  335644015  |  36	  |     9.8     | 
CRAM v3      |  308745955  |  24	  |     9.3     | 
CRAM v3+bz2  |  289888505  |  32 	  |     -       | Via Scramble -j
CRAM v3+lzma |  282989638  | 105 	  |     -       | Via Scramble -Z
CRAM v3 MAX  |  281666960  | 166 	  |     -       | Via Scramble -9 -jZ (bzip2, lzma)

Scramble was used to test bzip2, lzma and both combined along with
compression level 9 for maximum shrinkage.

