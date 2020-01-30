---
permalink: /howtos/headers.html
layout: default
title: Header viewing options
---
## SAMtools & BCFtools header viewing options

The roles of the `-h` and `-H` options in `samtools view` and `bcftools view` have historically
been inconsistent and confusing.
The 1.15 releases improve this by adding new `head` commands alongside the previous releases'
consistent sets of `view` long options.

To display only the headers of a SAM/BAM/CRAM or VCF/BCF file, use `head`:

```
samtools head FILE
bcftools head FILE
```

Both `head` commands also have options to limit the number of header lines displayed and/or to display the
first few alignment/variant records as well.

When displaying all or a subset of alignment/variant records, the two `view` commands differ as to
whether header lines will be included by default:

* `samtools view` omits headers by default; `-h` can be used to include them.

* `bcftools view` includes headers by default; `-H` can be used to omit them.

Both `view` commands accept the same long options specifying whether to include header lines:

```
samtools view [--no-header]   FILE		samtools view  --with-header|-h FILE
bcftools view  --no-header|-H FILE		bcftools view [--with-header]   FILE
```

The square brackets indicate options that are strictly speaking unnecessary as they select the
default behaviour anyway. However especially in scripts it can be useful to use them for explicitness,
and they have the effect of resetting previously specified `--header-only`/`--no-header`/`--with-header` options.

The `view` commands also have an option to display only headers, similarly to `head` above:

```
samtools view --header-only FILE
bcftools view --header-only FILE
```

For compatibility with earlier versions, there are also equivalent `view` short options.
On the command line we recommend using the more succinct `head` commands instead; trying to remember the
`view --header-only` short options mostly results only in confusion about which of the `-h`/`-H` options
does what with which command.
(Note that to be exactly equivalent to `head`, these `view` commands would also need to use `--no-PG`
or `--no-version` respectively.)

### Compatibility details

* For SAMtools: the `head` command is new in release 1.15, and the header-related `view` long options were new in 1.13.
* For BCFtools: the `head` command is new in release 1.15, and the `view --with-header` option was new in 1.14.
  (BCFtools view has had `--no-header` and `--header-only` since release 1.0.)
