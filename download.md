---
layout: default
title: SAMtools/BCFtools/HTSlib - Downloads
highlighting: yes
---
Current releases
----------------

**SAMtools** and **BCFtools** are distributed as individual packages.
The code uses HTSlib internally, but these source packages contain their own
copies of htslib so they can be built independently.

**HTSlib** is also distributed as a separate package which can be installed
if you are writing your own programs against the HTSlib API.
HTSlib also provides the **bgzip**, **htsfile**, and **tabix** utilities,
so you may also want to build and install HTSlib to get these utilities,
or see the additional instructions in `INSTALL` to install them from a
samtools or bcftools source package.

Download current source releases:
&nbsp;
<a class="btn btn-success" href="https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2"><span class="glyphicon glyphicon-download-alt" aria-hidden="true"></span> samtools-1.4</a>
&emsp;
<a class="btn btn-success" href="https://github.com/samtools/bcftools/releases/download/1.4/bcftools-1.4.tar.bz2"><span class="glyphicon glyphicon-download-alt" aria-hidden="true"></span> bcftools-1.4</a>
&emsp;
<a class="btn btn-success" href="https://github.com/samtools/htslib/releases/download/1.4/htslib-1.4.tar.bz2"><span class="glyphicon glyphicon-download-alt" aria-hidden="true"></span> htslib-1.4</a>

See also release notes for
[**samtools**](https://github.com/samtools/samtools/releases/latest/),
[**bcftools**](https://github.com/samtools/bcftools/releases/latest/),
and [**htslib**](https://github.com/samtools/htslib/releases/latest/).

New releases are announced on the [samtools mailing lists] and by [@htslib]
on Twitter.
Previous releases are available from the
[samtools GitHub organisation](https://github.com/samtools/)
(see [samtools](https://github.com/samtools/samtools/releases/),
[bcftools](https://github.com/samtools/bcftools/releases/),
or [htslib](https://github.com/samtools/htslib/releases/) releases)
or from the
[samtools Sourceforge project](http://sourceforge.net/projects/samtools/files/samtools/).

[@htslib]: https://twitter.com/htslib
[samtools mailing lists]: ../support#mailing-lists

### Building and installing

Building each desired package from source is very simple:

{% highlight sh %}
cd samtools-1.x    # and similarly for bcftools and htslib
make
make prefix=/where/to/install install
{% endhighlight %}

See `INSTALL` in each of the source directories for further details.

The executable programs will be installed to a `bin` subdirectory under
your specified prefix, so you may wish to add this directory to your $PATH:

{% highlight sh %}
export PATH=/where/to/install/bin:$PATH    # for sh or bash users
{% endhighlight %}
{% highlight csh %}
setenv PATH /where/to/install/bin:$PATH    # for csh users
{% endhighlight %}


Historical SAMtools/BCFtools 0.1.x releases
-------------------------------------------

Prior to the introduction of HTSlib, SAMtools and BCFtools were distributed
in a single samtools-0.1.x package.
These old versions remain available from the [Sourceforge samtools project](http://sourceforge.net/projects/samtools/files/samtools/).
