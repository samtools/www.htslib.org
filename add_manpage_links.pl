#!/usr/bin/perl

use strict;
use warnings;

my ($name) = @ARGV;
open(my $in, '<', $name) || die "Couldn't open $name : $!\n";
open(my $out, '>', "$name~") || die "Couldn't open $name~ for writing : $!\n";
while (<$in>) {
    if ($. == 1 && /^---$/) {
        print $out $_;
        add_frontmatter_redirects($in, $out);
    } elsif (/^<h1 id="SYNOPSIS">/) {
        print $out $_;
        add_samtools_link_synopsis($in, $out);
    } elsif (/^<h1 id="SEE_ALSO">/) {
        print $out $_;
        add_link_see_also($in, $out);
    } else {
        print $out $_;
    }
}
close($out) || die "Error closing $name~ : $!\n";
close($in) || die "Error closing $name : $!\n";
rename("$name~", $name) || die "Couldn't mv $name~ $name : $!\n";

sub add_frontmatter_redirects {
    my ($in, $out) = @_;

    my %redirects = (
        'bgzip.html'            => 'bgzip.1.html',
        'faidx.html'            => 'faidx.5.html',
        'htsfile.html'          => 'htsfile.1.html',
        'htslib-s3-plugin.html' => 'htslib-s3-plugin.7.html',
        'sam.html'              => 'sam.5.html',
        'samtools.html'         => 'samtools.1.html',
        'tabix.html'            => 'tabix.1.html',
        'vcf.html'              => 'vcf.5.html');

    my $redirect;
    while (<$in>) {
        if (m#^permalink:\s+(/doc(?:/[0-9][^/]+))/(\S+)# && exists($redirects{$2})) {
            $redirect = "${1}/$redirects{$2}";
        } elsif (m#^permalink:\s+(/doc(?:/[0-9][^/]+))/(samtools-[a-z]+)\.html#) {
            $redirect = "${1}/${2}.1.html";
        } elsif (/^---$/ && $redirect) {
            print $out "redirect_from: $redirect\n";
        }
        print $out $_;
        last if (/^---$/);
    }
}

sub add_samtools_link_synopsis {
    my ($in, $out) = @_;

    while (my $l = <$in>) {
        if ($l =~ /^samtools\s+([a-z]+)(.*)/) {
            my $extra = $2 || '';
            print $out qq[samtools <a href="samtools-$1.html">$1</a>$extra\n];
        } elsif ($l =~ /^<h1/) {
            print $out $l;
            last;
        } else {
            print $out $l;
        }
    }
}

sub add_link_see_also {
    my ($in, $out) = @_;

    my $page_names = qr/
        bcftools
        | bgzip
        | fadix
        | htsfile
        | htslib-s3-plugin
        | sam
        | samtools-?[a-z]*
        | tabix
        | vcf
        /x;

    while (my $l = <$in>) {
        if ($l =~ /^<em>($page_names)<\/em>\s*\((\d)\)(,?)/) {
            my $comma = $3 || '';
            print $out qq[<a href="$1.html"><em>$1</em></a> ($2)$comma\n];
        } elsif ($l =~ /^<em>[^<]+<\/em>/) {
            print $out $l;
        } else {
            print $out $l;
            last;
        }
    }
}
