---
layout: default
title: Zlib Benchmarks
highlighting: yes
---
A comparison of Zlib implementations
------------------------------------

Samtools has been minimally tested against both the Intel-optimised and
CloudFlare-optimised zlibs and shown to work.

They can be downloaded from:

* <https://github.com/jtkukunas/zlib>     # Intel
* <https://github.com/cloudflare/zlib>    # CloudFlare


External benchmarks comparing the various zlibs (on non-BAM data) are
available at: 

* <https://www.snellman.net/blog/archive/2014-08-04-comparison-of-intel-and-cloudflare-zlib-patches.html>

It is recommended that you perform your own rigorous tests for an
entire pipeline if you wish to switch to one of the optimised zlib
implementations, but we present some basic read and read/write tests
of our own.

In addition to these two full Zlib compatible implementations there is
an Intel igzip library.  This has a different and incompatible API,
but offers a better tradeoff between size and time than either of
these libraries.  It needs more work to make it seemlessly fit in to
Htslib and Samtools, but a branch containing this is benchmarked to
demonstrate the gains and losses.

### Reading (flagstats)

(1 thread, elapsed ~= cpu time)

{:.table}
Zlib Implementation | CPU Time
------------------- | :------:
Original            | 1m 17s
Intel               | 1m 17s
CloudFlare          | 1m 17s

### Re-writing (BAM->BAM)

(4 threads, fastest elapsed/CPU, size)

{:.table}
Zlib Implementation | Level | Elapsed Time  | CPU Time  |       Size
------------------- | :---: | :-----------: | --------: | :--------:
Original            |     6 | 4m 31s        | 12m 50s   | 6499959405
Intel               |     6 | 3m 49s        |  9m 57s   | 6580893700
CloudFlare          |     6 | 3m 24s        |  8m 13s   | 6445044187

At minimum compression level:

{:.table}
Zlib Implementation | Level | Elapsed Time  | CPU Time  |       Size
------------------- | :---: | :-----------: | --------: | :--------:
Original            |     1 | 3m 00s        |  6m 42s   | 7042519949
Intel               |     1 | 2m 24s        |  4m 15s   | 9019497823
CloudFlare          |     1 | 2m 53s        |  6m 12s   | 6816769584
igzip               |     1 | 2m 09s        |  3m 26s   | 7324800394


### Summary


CloudFlare zlib is faster and smaller than the default implementation
at both standard compression level and also level -1.

The Intel zlib specialises in fastest compression at level -1, albeit
at a weaker compression ratio.

Both can be used with LD_PRELOAD or LD_LIBRARY_PATH without needing to
recompile samtools/htslib. (Or samtools could be built linking against
one of these libraries and adding a -rpath, but that needs manual
modification of the Makefile.)

The igzip specialised BAM compression code runs fastest for zlib-1
equivalent while giving an acceptable compression ratio, but code
integration is currently not simple.  Arguably using a different
format entirely, such as LZ4, would be preferable for temporary
intermediate files. (Not benchmarked here.)


Full reading benchmarks
-----------------------

### Original ZLib

{% highlight sh %}
$ time ./samtools flagstat ~/scratch/data/9827_2#49.bam
56463236 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
269248 + 0 duplicates
55357963 + 0 mapped (98.04% : N/A)
56463236 + 0 paired in sequencing
28231618 + 0 read1
28231618 + 0 read2
54363468 + 0 properly paired (96.28% : N/A)
55062652 + 0 with itself and mate mapped
295311 + 0 singletons (0.52% : N/A)
360264 + 0 with mate mapped to a different chr
300699 + 0 with mate mapped to a different chr (mapQ>=5)

real    1m16.934s
user    1m15.520s
sys     0m1.300s

$ time ./samtools flagstat ~/scratch/data/9827_2#49.bam
56463236 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
269248 + 0 duplicates
55357963 + 0 mapped (98.04% : N/A)
56463236 + 0 paired in sequencing
28231618 + 0 read1
28231618 + 0 read2
54363468 + 0 properly paired (96.28% : N/A)
55062652 + 0 with itself and mate mapped
295311 + 0 singletons (0.52% : N/A)
360264 + 0 with mate mapped to a different chr
300699 + 0 with mate mapped to a different chr (mapQ>=5)

real    1m16.891s
user    1m15.820s
sys     0m0.950s
{% endhighlight %}

### Intel ZLib

{% highlight sh %}
$ time LD_PRELOAD=/tmp/zlib.intel/libz.so ./samtools flagstat ~/scratch/data/9827_2#49.bam
56463236 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
269248 + 0 duplicates
55357963 + 0 mapped (98.04% : N/A)
56463236 + 0 paired in sequencing
28231618 + 0 read1
28231618 + 0 read2
54363468 + 0 properly paired (96.28% : N/A)
55062652 + 0 with itself and mate mapped
295311 + 0 singletons (0.52% : N/A)
360264 + 0 with mate mapped to a different chr
300699 + 0 with mate mapped to a different chr (mapQ>=5)

real    1m17.579s
user    1m16.140s
sys     0m1.310s

$ time LD_PRELOAD=/tmp/zlib.intel/libz.so ./samtools flagstat ~/scratch/data/9827_2#49.bam
56463236 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
269248 + 0 duplicates
55357963 + 0 mapped (98.04% : N/A)
56463236 + 0 paired in sequencing
28231618 + 0 read1
28231618 + 0 read2
54363468 + 0 properly paired (96.28% : N/A)
55062652 + 0 with itself and mate mapped
295311 + 0 singletons (0.52% : N/A)
360264 + 0 with mate mapped to a different chr
300699 + 0 with mate mapped to a different chr (mapQ>=5)

real    1m17.280s
user    1m15.960s
sys     0m1.200s
{% endhighlight %}

### CloudFlare ZLib

{% highlight sh %}
$ time LD_PRELOAD=/tmp/zlib.cloudflare/libz.so ./samtools flagstat ~/scratch/data/9827_2#49.bam
56463236 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
269248 + 0 duplicates
55357963 + 0 mapped (98.04% : N/A)
56463236 + 0 paired in sequencing
28231618 + 0 read1
28231618 + 0 read2
54363468 + 0 properly paired (96.28% : N/A)
55062652 + 0 with itself and mate mapped
295311 + 0 singletons (0.52% : N/A)
360264 + 0 with mate mapped to a different chr
300699 + 0 with mate mapped to a different chr (mapQ>=5)

real    1m17.375s
user    1m16.310s
sys     0m0.940s

$ time LD_PRELOAD=/tmp/zlib.cloudflare/libz.so ./samtools flagstat ~/scratch/data/9827_2#49.bam
56463236 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
269248 + 0 duplicates
55357963 + 0 mapped (98.04% : N/A)
56463236 + 0 paired in sequencing
28231618 + 0 read1
28231618 + 0 read2
54363468 + 0 properly paired (96.28% : N/A)
55062652 + 0 with itself and mate mapped
295311 + 0 singletons (0.52% : N/A)
360264 + 0 with mate mapped to a different chr
300699 + 0 with mate mapped to a different chr (mapQ>=5)

real    1m17.144s
user    1m16.020s
sys     0m1.000s
{% endhighlight %}

Full read/write benchmarks
--------------------------

### Original ZLib, default compression level

{% highlight sh %}
$ time ./samtools view -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    4m33.612s
user    12m38.000s
sys     0m19.920s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 6499959405 Apr 27 11:42 /tmp/_.bam

$ time ./samtools view -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    4m31.879s
user    12m29.380s
sys     0m20.930s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 6499959405 Apr 27 11:53 /tmp/_.bam
{% endhighlight %}

### Intel ZLib, default compression level

{% highlight sh %}
$ time LD_PRELOAD=/tmp/zlib.intel/libz.so ./samtools view -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    3m49.306s
user    9m36.150s
sys     0m21.290s
$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 6580893700 Apr 27 11:45 /tmp/_.bam

$ time LD_PRELOAD=/tmp/zlib.intel/libz.so ./samtools view -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    3m52.742s
user    9m43.360s
sys     0m20.590s
$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 6580893700 Apr 27 11:57 /tmp/_.bam
{% endhighlight %}

### CloudFlare ZLib, default compression level

{% highlight sh %}
$ time LD_PRELOAD=/tmp/zlib.cloudflare/libz.so ./samtools view -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    3m23.885s
user    7m51.620s
sys     0m21.210s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 6445044187 Apr 27 11:49 /tmp/_.bam

$ time LD_PRELOAD=/tmp/zlib.cloudflare/libz.so ./samtools view -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    3m24.938s
user    7m55.610s
sys     0m21.430s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 6445044187 Apr 27 12:01 /tmp/_.bam
{% endhighlight %}

### Original ZLib, lowest compression level

{% highlight sh %}
$ time ./samtools view -1 -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    2m59.867s
user    6m21.680s
sys     0m20.210s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 7042519949 Apr 27 11:23 /tmp/_.bam

$ time ./samtools view -1 -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    3m1.529s
user    6m21.060s
sys     0m22.180s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 7042519949 Apr 27 11:32 /tmp/_.bam
{% endhighlight %}

### Intel ZLib, lowest compression level

{% highlight sh %}
$ time LD_PRELOAD=/tmp/zlib.intel/libz.so ./samtools view -1 -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    2m24.214s
user    4m0.170s
sys     0m14.790s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 9019497823 Apr 27 11:26 /tmp/_.bam

$ time LD_PRELOAD=/tmp/zlib.intel/libz.so ./samtools view -1 -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    2m25.409s
user    4m2.690s
sys     0m13.300s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 9019497823 Apr 27 11:34 /tmp/_.bam
{% endhighlight %}

### CloudFlare ZLib, lowest compression level

{% highlight sh %}
$ time LD_PRELOAD=/tmp/zlib.cloudflare/libz.so ./samtools view -1 -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    2m52.685s
user    5m39.080s
sys     0m23.250s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 6816769584 Apr 27 11:29 /tmp/_.bam

$ time LD_PRELOAD=/tmp/zlib.cloudflare/libz.so ./samtools view -1 -b -@ 4 -o /tmp/_.bam /nfs/users/nfs_j/jkb/scratch/data/9827_2#49.bam

real    2m54.839s
user    5m44.820s
sys     0m23.270s

$ ls -l /tmp/_.bam
-rw-r--r-- 1 jkb team117 6816769584 Apr 27 11:37 /tmp/_.bam
{% endhighlight %}

### Intel igzip library, lowest compresion level

{% highlight sh %}
$ time ./test/test_view -@ 4 -l 1 -b ~/scratch/data/9827_2#49.bam > /tmp/_.bam;ls -l /tmp/_.bam

real    2m9.550s
user    3m17.260s
sys     0m9.190s
-rw-r--r-- 1 jkb team117 7324800394 Apr 28 14:54 /tmp/_.bam

$ time ./test/test_view -@ 4 -l 1 -b ~/scratch/data/9827_2#49.bam > /tmp/_.bam;ls -l /tmp/_.bam

real    2m14.863s
user    3m17.730s
sys     0m7.990s
-rw-r--r-- 1 jkb team117 7324800394 Apr 28 14:59 /tmp/_.bam
{% endhighlight %}
