HTSLIB   = ../htslib
SAMTOOLS = ../samtools
BCFTOOLS = ../bcftools

ifeq ($(VERSION),)
HTSLIB_VERSION ?= $(shell git --git-dir=$(HTSLIB)/.git describe --match '[0-9].[0-9]*' | sed 's/-.*//')
SAMTOOLS_VERSION ?= $(shell git --git-dir=$(SAMTOOLS)/.git describe --match '[0-9].[0-9]*' | sed 's/-.*//')
BCFTOOLS_VERSION ?= $(shell git --git-dir=$(BCFTOOLS)/.git describe --match '[0-9].[0-9]*' | sed 's/-.*//')
else
HTSLIB_VERSION ?= $(VERSION)
SAMTOOLS_VERSION ?= $(VERSION)
BCFTOOLS_VERSION ?= $(VERSION)
COMMIT ?= $(VERSION)
endif

ifeq ($(COMMIT),)
HTSLIB_COMMIT ?= $(HTSLIB_VERSION)
SAMTOOLS_COMMIT ?= $(SAMTOOLS_VERSION)
BCFTOOLS_COMMIT ?= $(BCFTOOLS_VERSION)
else
HTSLIB_COMMIT ?= $(COMMIT)
SAMTOOLS_COMMIT ?= $(COMMIT)
BCFTOOLS_COMMIT ?= $(COMMIT)
endif

all:	doc.md doc/samtools.html doc/bcftools.html \
	doc/htsfile.html doc/tabix.html \
	doc/faidx.html doc/sam.html doc/vcf.html

doc.md: doc/samtools.html doc/bcftools.html doc/htsfile.html \
	doc/tabix.html doc/faidx.html doc/sam.html doc/vcf.html
	./update_doc_md.pl

doc/samtools.html: doc/samtools-$(SAMTOOLS_VERSION).html
	sed '/^permalink:/s/-$(SAMTOOLS_VERSION)//' $< > $@

doc/bcftools.html: doc/bcftools-$(BCFTOOLS_VERSION).html
	sed '/^permalink:/s/-$(BCFTOOLS_VERSION)//' $< > $@

doc/htsfile.html: doc/htsfile-$(HTSLIB_VERSION).html
	sed '/^permalink:/s/-$(HTSLIB_VERSION)//' $< > $@

doc/tabix.html: doc/tabix-$(HTSLIB_VERSION).html
	sed '/^permalink:/s/-$(HTSLIB_VERSION)//' $< > $@

GRIND = man2fhtml --mode jekyll --location /$@ --output $@

doc/samtools-$(SAMTOOLS_VERSION).html:
	git --git-dir=$(SAMTOOLS)/.git show $(SAMTOOLS_COMMIT):samtools.1 | $(GRIND)

BCFTOOLS_DOC_DATE = $(shell git --git-dir=$(BCFTOOLS)/.git log -n 1 0.1.0..$(BCFTOOLS_COMMIT) --date=short --pretty=format:%cd -- doc/bcftools.txt)

doc/bcftools-$(BCFTOOLS_VERSION).html:
	git --git-dir=$(BCFTOOLS)/.git show $(BCFTOOLS_COMMIT):doc/bcftools.txt > doc/bcftools-$(BCFTOOLS_VERSION).txt && \
	a2x  -adate='$(BCFTOOLS_DOC_DATE)' -aversion=$(BCFTOOLS_VERSION) --doctype manpage --format xhtml doc/bcftools-$(BCFTOOLS_VERSION).txt && \
	rm doc/bcftools-$(BCFTOOLS_VERSION).txt

doc/htsfile-$(HTSLIB_VERSION).html:
	git --git-dir=$(HTSLIB)/.git show $(HTSLIB_COMMIT):htsfile.1 | $(GRIND)

doc/tabix-$(HTSLIB_VERSION).html:
	git --git-dir=$(HTSLIB)/.git show $(HTSLIB_COMMIT):tabix.1 | $(GRIND)

doc/faidx.html: $(HTSLIB)/faidx.5
	$(GRIND) $<

doc/sam.html: $(HTSLIB)/sam.5
	$(GRIND) $<

doc/vcf.html: $(HTSLIB)/vcf.5
	$(GRIND) $<
