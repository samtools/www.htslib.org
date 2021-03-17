#!/usr/bin/perl

my @in = @ARGV;
@in = sort by_version @in;
print "@in";

sub by_version {
    my @as = split(/\./, $a);
    my @bs = split(/\./, $b);

    while (@as && @bs) {
	my $an = shift(@as);
	my $bn = shift(@bs);
	if ($an =~ /^\d+$/ && $bn =~ /^\d+$/) {
	    return $an <=> $bn if ($an <=> $bn);
	} else {
	    return $an cmp $bn if ($an cmp $bn);
	}
    }
    if (!@as) { return -1; }
    if (!@bs) { return  1; }
    return 0;
}
