#!/usr/bin/perl
use strict;
use warnings;

sub dotted_version {
    my ($vers) = @_;

    die "vers = $vers" unless (length($vers) == 9);
    if (substr($vers, 6, 3) eq '000') {
	sprintf("%d.%d", substr($vers,0, 3), substr($vers, 3, 3));
    } else {
	return sprintf("%d.%d.%d", substr($vers,0, 3),
		       substr($vers, 3, 3), substr($vers, 6, 3));
    }
}

my $doc_dir = 'doc';
my $doc_md = 'doc.md';

# Scan doc directory for versioned man pages
my %manpages;
opendir(my $dir, $doc_dir) || die "Couldn't open $doc_dir : $!\n";
while ($_ = readdir($dir)) {
    if (/^(bcftools|htsfile|samtools|tabix|bgzip)-(\d+)\.(\d+)(?:.(\d+))?\.html/) {
	my $v = sprintf("%03d%03d%03d", $2, $3, $4 ? $4 : 0);
	$manpages{$1}->{$v} = $_;
	if ($1 eq 'samtools' && $v eq '000001019') {
	    # Pre-htslib, samtools and bcftools shared a repository
	    $manpages{bcftools}->{$v} = $_;
	}
    }
}
closedir($dir) || die "Error reading $doc_dir : $!\n";

# bgzip and tabix shared a page before release 1.8
foreach my $v (keys %{$manpages{tabix}}) {
    if ($v le '001008000' && !exists($manpages{bgzip}->{$v})) {
	$manpages{bgzip}->{$v} = $manpages{tabix}->{$v};
    }
}

# Rewrite doc.md file with new man page links
open(my $doc_in, '<', $doc_md) || die "Couldn't open $doc_md : $!\n";
open(my $doc_out, '>', "$doc_md~")
    || die "Couldn't open $doc_md~ for writing : $!\n";
my $skip = 0;
my $manpages_section = 0;
while (<$doc_in>) {
    if (/^##\s+/) {
	$manpages_section = /Manual pages/;
	if ($skip && !$manpages_section) { $skip = 0; }
    }
    next if ($skip);
    if (!$manpages_section || !/^\* \[/) {
	print $doc_out $_;
	next;
    }

    # This is the section we want to replace
    $skip = 1;

    foreach my $tool (qw(bcftools bgzip htsfile samtools tabix)) {
	my $page = $tool;

	# Get versions, newest to oldest
	my @versions = sort { $b cmp $a } keys %{$manpages{$page}};

	# For newest, we use the unversioned page.
	my $newest = shift(@versions);
	my $v = dotted_version($newest);
	if (@versions) {
	    print $doc_out "* [$tool $v]($page.html) (older versions:\n";
	} else {
	    print $doc_out "* [$tool $v]($page.html)\n";
	}

	# Add list of old manpages
	while (my $vers = shift @versions) {
	    my $v = dotted_version($vers);
	    my $note = (($tool eq 'bcftools' && $v eq '0.1.19')
			? ' "included in samtools-0.1.19"' : '');
	    my $trail = @versions ? ',' : ')';
	    print $doc_out "      [$v]($manpages{$page}->{$vers}$note)$trail\n";
	}
    }
    print $doc_out "\n";
}
close($doc_in) || die "Error reading $doc_md : $!\n";
close($doc_out) || die "Error writing to $doc_md~ : $!\n";

rename("$doc_md~", "$doc_md")
    || die "Couldn't rename $doc_md~ to $doc_md : $!\n";
