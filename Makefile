HTSLIB   = ../htslib
SAMTOOLS = ../samtools
BCFTOOLS = ../bcftools

COMMIT ?= $(VERSION)

all:	doc/samtools-$(VERSION).html \
	doc/bcftools-$(VERSION).html \
	doc/htsfile-$(VERSION).html \
	doc/tabix-$(VERSION).html \
	doc/faidx.html doc/sam.html doc/vcf.html

GRIND = man2fhtml --mode jekyll --location /$@ --output $@

doc/samtools-$(VERSION).html:
	git --git-dir=$(SAMTOOLS)/.git show $(COMMIT):samtools.1 | $(GRIND)

doc/bcftools-$(VERSION).html:
	git --git-dir=$(BCFTOOLS)/.git show $(COMMIT):doc/bcftools.html > $@

doc/htsfile-$(VERSION).html:
	git --git-dir=$(HTSLIB)/.git show $(COMMIT):htsfile.1 | $(GRIND)

doc/tabix-$(VERSION).html:
	git --git-dir=$(HTSLIB)/.git show $(COMMIT):tabix.1 | $(GRIND)

doc/faidx.html: $(HTSLIB)/faidx.5
	$(GRIND) $<

doc/sam.html: $(HTSLIB)/sam.5
	$(GRIND) $<

doc/vcf.html: $(HTSLIB)/vcf.5
	$(GRIND) $<
