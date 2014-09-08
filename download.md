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
HTSlib also provides the **bgzip** and **tabix** utilities, so you may
also want to build and install HTSlib to get these utilities, or see the
additional instructions in `INSTALL` to install them from a samtools or
bcftools source package.

The current release of each package is **1.0**.

<a class="btn btn-primary" href="http://sourceforge.net/projects/samtools/files/samtools/1.0/"><i class="glyphicon glyphicon-save"></i> Download source releases here</a>

See release notes for [bcftools]({{ site.baseurl }}/bcftools_release_notes).
<!-- TODO Make releases available as GitHub releases in the three repos -->

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
