#!/usr/bin/perl -w

use strict;
use warnings;

die("Usage: src2html.pl <in1.c> [<in2.c> [...]]\n") if (@ARGV == 0);

for my $fn (@ARGV) {
  my $alt_fn = $fn;
  $alt_fn =~ s/.*\/([^\s\/]+)$/$1/;
  my ($fhout, $fhin);
  open($fhout, ">$fn.html") || die;
  print $fhout qq(<html>
<head><title>$alt_fn</title></head>
<body>
<textarea name="code" class="c" cols="60" rows="10">
);

  open($fhin, $fn) || die;
  while (<$fhin>) {
	print $fhout $_;
  }
  close($fhin);

print $fhout qq(
</textarea>
<link type="text/css" rel="stylesheet" href="../SyntaxHighlighter.css"></link>
<script language="javascript" src="../shCore.js"></script>
<script language="javascript" src="../shBrushCpp.js"></script>
<script language="javascript">
  dp.SyntaxHighlighter.HighlightAll('code');
</script>
</body>
</html>
);
  close($fhout);
}
