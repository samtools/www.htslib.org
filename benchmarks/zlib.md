---
permalink: /benchmarks/zlib.html
layout: default
title: Zlib Benchmarks
highlighting: yes
---
A comparison of Zlib implementations
------------------------------------

(Updated Sep 2019)

<!-- Tested on my deskpro (Intel(R) Core(TM) i5-4570 CPU @ 3.20GHz)
using ~/scratch/data/9827_2#49.bam:

READ: perf stat ./test/test_view -B ~/scratch/data/9827_2#49.bam
R/W:  perf stat ./test/test_view -@4 -b ~/scratch/data/9827_2#49.bam |wc -c
-->

Htslib / Samtools have been tested against the vanilla Zlib (1.2.11)
and the the Intel-optimised and CloudFlare-optimised zlibs and shown
to work.  Htslib also supports, indeed recommends, using libdeflate.
The benchmarks below indicate why this is recommended.

They can be downloaded from:

* <https://github.com/jtkukunas/zlib>      # Intel
* <https://github.com/cloudflare/zlib>     # CloudFlare
* <https://github.com/ebiggers/libdeflate> # Libdeflate

External benchmarks comparing the various zlibs (on non-BAM data) are
available at:

* <https://www.snellman.net/blog/archive/2014-08-04-comparison-of-intel-and-cloudflare-zlib-patches.html>

Previously we also tested Intel's igzip library on an older system.
Like libdeflate this has a different and incompatible API to standard
Zlib.  Benchmarks for this are estimated by applying a scaling
(reduction) factor to the old benchmarks based on observed speedup for
the other libraries.

### Reading (htslib's test_view -B in.bam)

(1 thread, elapsed ~= cpu time)

{:.table}
Zlib Implementation | CPU Time
------------------- | :------:
Original            | 1m 17s
Intel               | 1m 15s
CloudFlare          | 1m  4s
Libdeflate          | 0m 43s

### Re-writing (BAM->BAM; test_view -b -@4 in.bam | wc -c)

(4 threads, fastest elapsed/CPU, size)

{:.table}
Zlib Implementation | Level | Elapsed Time  | CPU Time  |       Size
------------------- | :---: | :-----------: | --------: | :--------:
Original            |     6 | 3m 12s        | 11m  6s   | 6499959405
Intel               |     6 | 2m 11s        |  8m 26s   | 6580893700
CloudFlare          |     6 | 1m 45s        |  6m 38s   | 6445044187
Libdeflate	    |     6 | 1m 20s        |  4m 57s   | 6526400042

At minimum compression level:

{:.table}
Zlib Implementation | Level | Elapsed Time  | CPU Time  |       Size
------------------- | :---: | :-----------: | --------: | :--------:
Original            |     1 | 1m 31s        |  5m 43s   | 7042519949
Intel               |     1 | 0m 55s        |  3m 21s   | 9019497823
CloudFlare          |     1 | 1m 14s        |  4m 40s   | 6816769584
Libdeflate          |     1 | 0m 58s        |  3m 26s   | 6863351382
igzip (est.)        |     1 | 0m 50s        |  2m 39s   | 7324800394


### Summary

All implementations fair favourably compared to the vanilla zlib
implementation.

Libdeflate is our recommended library as can be seen above.  It offers
a considerable speed improvement for both decoding and encoding. File
sizes are a bit larger than the default zlib, although for
(considerably) more CPU time it scales up to level 12 compared the
default level 9 which offers the best size (not shown here).

CloudFlare's zlib offers a good default file size for not much more
time than libdeflate so is still a worthy contender.

The Intel zlib specialises in fast compression at level -1, albeit
at a weaker compression ratio.  The BAM-specific mode of igzip takes
this further offering the best speeds, but file sizes are still behind
most other implementations.  It has potential use for intermediate
temporary files, like likely something completely different such as
LZ4 would be preferable.

Both Intel and CloudFlare Zlib's can be used with LD_PRELOAD or
LD_LIBRARY_PATH without needing to recompile samtools/htslib, provided
it was built without enabling --with-libdeflate.
