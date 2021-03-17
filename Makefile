HTSLIB   = ../htslib
SAMTOOLS = ../samtools
BCFTOOLS = ../bcftools

MAN2FHTML = man2fhtml
ADD_MANPAGE_LINKS = ./add_manpage_links.pl
VERS_SORT = ./vers_sort.pl

all:	samtools-doc htslib-doc bcftools-doc update_doc.md

samtools-doc:
	@ for i in $(SAMTOOLS)/doc/*.1; do \
	    base=`echo $$i | sed 's:.*/::;s:\.[1-9]$$::'`; \
	    echo Processing $$i;\
	    $(MAN2FHTML) --mode jekyll --location /doc/$$base.html --output doc/$$base.html < $$i;\
	    $(ADD_MANPAGE_LINKS) doc/$$base.html;\
	done

BCFTOOLS_VERSION ?= $(shell git --git-dir=$(BCFTOOLS)/.git describe --match '[0-9].[0-9]*' | sed 's/-.*//')
BCFTOOLS_DOC_DATE = $(shell git --git-dir=$(BCFTOOLS)/.git log -n 1 0.1.0.. --date=short --pretty=format:%cd -- doc/bcftools.txt)

bcftools-doc:
	a2x -adate='$(BCFTOOLS_DOC_DATE)' -aversion=$(BCFTOOLS_VERSION) --doctype manpage --format xhtml -D doc $(BCFTOOLS)/doc/bcftools.txt

htslib-doc:
	@ for i in $(HTSLIB)/*.[1-9]; do \
	    case $$i in \
	    *".so."*) \
	        ;; \
	    *) \
	        base=`echo $$i | sed 's:.*/::;s:\.[1-9]$$::'`; \
	        echo Processing $$i;\
	        $(MAN2FHTML) --mode jekyll --location /doc/$$base.html --output doc/$$base.html < $$i;\
	        $(ADD_MANPAGE_LINKS) doc/$$base.html;\
	        ;; \
	    esac \
	done

update_doc.md:
	vers="";for v in `$(VERS_SORT) doc/[0-9]*|sed 's#doc/##g'`;do vers="$$vers$${vers:+, }[$$v]($$v)";done; \
	mv doc.md doc.md~ && sed "s#for releases:.*#for releases: $$vers#" doc.md~ > doc.md
	rm doc.md~

