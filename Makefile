HTSLIB   = ../htslib
SAMTOOLS = ../samtools
BCFTOOLS = ../bcftools

MAN2FHTML = ./man2fhtml

all:	samtools-doc htslib-doc bcftools-doc update_doc.md

# Note this may need running twice to get SYNOPSIS links correct
# if we removed doc/*.html first. (So don't do that.)
samtools-doc:
	@ for i in $(SAMTOOLS)/doc/*.1; do \
	    base=`echo $$i | sed 's:.*/::'`; \
	    echo Processing $$i;\
	    $(MAN2FHTML) --mode jekyll --location /doc/$$base.html --output doc/$$base.html < $$i;\
	done

BCFTOOLS_VERSION ?= $(shell git --git-dir=$(BCFTOOLS)/.git describe --match '[0-9].[0-9]*' | sed 's/-.*//')
BCFTOOLS_DOC_DATE = $(shell git --git-dir=$(BCFTOOLS)/.git log -n 1 0.1.0.. --date=short --pretty=format:%cd -- doc/bcftools.txt)

bcftools-doc:
	a2x -adate='$(BCFTOOLS_DOC_DATE)' -aversion=$(BCFTOOLS_VERSION) --doctype manpage --format xhtml -D doc $(BCFTOOLS)/doc/bcftools.txt
	mv doc/bcftools.html doc/bcftools.1.html

htslib-doc:
	@ for i in $(HTSLIB)/*.[1-9]; do \
	    base=`echo $$i | sed 's:.*/::'`; \
	    echo Processing $$i;\
	    $(MAN2FHTML) --mode jekyll --location /doc/$$base.html --output doc/$$base.html < $$i;\
	done

update_doc.md:
	vers="";for v in `echo doc/[0-9]*|sed 's#doc/##g'`;do vers="$$vers$${vers:+, }[$$v]($$v)";done; \
	mv doc.md doc.md~ && sed "s#for releases:.*#for releases: $$vers#" doc.md~ > doc.md
	rm doc.md~

