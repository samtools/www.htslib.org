HTSLIB   = ../htslib
SAMTOOLS = ../samtools
BCFTOOLS = ../bcftools

MAN2FHTML = man2fhtml
ADD_MANPAGE_LINKS = ./add_manpage_links.pl
VERS_SORT = ./vers_sort.pl

ifeq ($(VERSION),)
HTSLIB_VERSION ?= $(shell git --git-dir=$(HTSLIB)/.git describe --match '[0-9].[0-9]*' | sed 's/-.*//')
SAMTOOLS_VERSION ?= $(shell git --git-dir=$(SAMTOOLS)/.git describe --match '[0-9].[0-9]*' | sed 's/-.*//')
BCFTOOLS_VERSION ?= $(shell git --git-dir=$(BCFTOOLS)/.git describe --match '[0-9].[0-9]*' | sed 's/-.*//')
else
HTSLIB_VERSION ?= $(VERSION)
SAMTOOLS_VERSION ?= $(VERSION)
BCFTOOLS_VERSION ?= $(VERSION)
endif

HTSLIB_VERSIONED_DOC ?= doc/$(HTSLIB_VERSION)
SAMTOOLS_VERSIONED_DOC ?= doc/$(SAMTOOLS_VERSION)
BCFTOOLS_VERSIONED_DOC ?= doc/$(BCFTOOLS_VERSION)

all:	samtools-doc htslib-doc bcftools-doc update_doc.md

samtools-doc: | $(SAMTOOLS_VERSIONED_DOC)
	@ for i in $(SAMTOOLS)/doc/*.1; do \
	    base=`echo $$i | sed 's:.*/::;s:\.[1-9]$$::'`; \
	    if [ "$$base" = "samtools-fastq" ] ; then continue ; fi ; \
	    echo Processing $$i;\
	    $(MAN2FHTML) --mode jekyll --location /$(SAMTOOLS_VERSIONED_DOC)/$$base.html --output $(SAMTOOLS_VERSIONED_DOC)/$$base.html < $$i;\
	    $(ADD_MANPAGE_LINKS) $(SAMTOOLS_VERSIONED_DOC)/$$base.html;\
	    sed -E '/^(permalink|redirect_from):/s#doc/$(SAMTOOLS_VERSION)/#doc/#' $(SAMTOOLS_VERSIONED_DOC)/$$base.html > doc/$$base.html ; \
	done && \
	cp doc/samtools-fastq.html $(SAMTOOLS_VERSIONED_DOC)/samtools-fastq.html

BCFTOOLS_DOC_DATE = $(shell git --git-dir=$(BCFTOOLS)/.git log -n 1 0.1.0.. --date=short --pretty=format:%cd -- doc/bcftools.txt)

bcftools-doc: | $(BCFTOOLS_VERSIONED_DOC)
	asciidoctor -adate='$(BCFTOOLS_DOC_DATE)' -aversion=$(BCFTOOLS_VERSION) --doctype manpage --backend xhtml -D doc $(BCFTOOLS)/doc/bcftools.txt && \
	sed 's#href="docbook-xsl\.css"#href="../docbook-xsl.css"#' doc/bcftools.html > $(BCFTOOLS_VERSIONED_DOC)/bcftools.html 

htslib-doc: | $(HTSLIB_VERSIONED_DOC)
	@ for i in $(HTSLIB)/*.[1-9]; do \
	    case $$i in \
	    *".so."*) \
	        ;; \
	    *) \
	        base=`echo $$i | sed 's:.*/::;s:\.[1-9]$$::'`; \
	        echo Processing $$i;\
	        $(MAN2FHTML) --mode jekyll --location /$(HTSLIB_VERSIONED_DOC)/$$base.html --output $(HTSLIB_VERSIONED_DOC)/$$base.html < $$i;\
	        $(ADD_MANPAGE_LINKS) $(HTSLIB_VERSIONED_DOC)/$$base.html;\
	        sed -E '/^(permalink|redirect_from):/s#doc/$(HTSLIB_VERSION)/#doc/#' $(HTSLIB_VERSIONED_DOC)/$$base.html > doc/$$base.html ; \
	        ;; \
	    esac \
	done

$(HTSLIB_VERSIONED_DOC) $(filter-out $(HTSLIB_VERSIONED_DOC),$(SAMTOOLS_VERSIONED_DOC)) $(filter-out $(HTSLIB_VERSIONED_DOC) $(SAMTOOLS_VERSIONED_DOC),$(BCFTOOLS_VERSIONED_DOC)):
	mkdir -p "$@" && \
	dir="$@" && \
	vers=$${dir##*/} && \
	printf -- '---\nlayout: default\ntitle: Samtools - Documentation\n---\n## Manual pages\n\nDocumentation for BCFtools, SAMtools, and HTSlib'"'"'s utilities is available\nby using <code>man <em>command</em></code> on the command line.\n' > $@/index.md && \
	printf 'The manual pages for the %s release are listed below.\n\n' "$$vers" >> $@/index.md && \
	( if [ "$(BCFTOOLS_VERSIONED_DOC)" = "$@" ] ; then \
	    printf '* [bcftools](bcftools.html)\n' >> "$@"/index.md ; \
	fi ; \
	if [ "$(HTSLIB_VERSIONED_DOC)" = "$@" ] ; then \
	    printf '* [bgzip](bgzip.html)\n* [htsfile](htsfile.html)\n' >> "$@"/index.md ; \
	fi ; \
	if [ "$(SAMTOOLS_VERSIONED_DOC)" = "$@" ] ; then \
	    printf '* [samtools](samtools.html)\n' >> "$@"/index.md ; \
	fi ; \
	if [ "$(HTSLIB_VERSIONED_DOC)" = "$@" ]	; then \
            printf '* [tabix](tabix.html)\n' >> "$@"/index.md ; \
        fi )

update_doc.md:
	vers=""; \
	for v in `$(VERS_SORT) doc/[0-9]*|sed 's#doc/##g'`; do \
	    if [ "doc/$$v" != "$(HTSLIB_VERSIONED_DOC)" ]   || \
	       [ "doc/$$v" != "$(SAMTOOLS_VERSIONED_DOC)" ] || \
	       [ "doc/$$v" != "$(BCFTOOLS_VERSIONED_DOC)" ] ; then \
	        vers="$$vers$${vers:+, }[$$v]($$v)"; \
	    fi ; \
	done; \
	mv doc.md doc.md~ && sed "s#for releases:.*#for releases: $$vers#" doc.md~ > doc.md
	rm doc.md~

