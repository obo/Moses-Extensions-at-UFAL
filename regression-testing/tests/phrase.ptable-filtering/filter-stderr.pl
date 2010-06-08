#!/usr/bin/perl
$x=0;
while (<>) {
  chomp;

	if (/^\[.* ; 2-2\]$/o) {
    my @lines;
    my $done = 0;
    while (!$done) {
      $x = <>;
      if ($x =~ /^\s*$/o) { $done = 1; } else {
        chomp $x;
        $x =~ s/^\s+//o;
        push @lines, $x;
      }
    }
    my $c = 0;
    foreach my $x (sort @lines) {
      $c++;
      print "TRANSLATION_OPTION_$c=$x\n";
    }
  }

  next unless /^BEST TRANSLATION:/;
  s/^BEST TRANSLATION:\s*//;
  s/\s*\[111+.*$//;
  $x++;
  print "TRANSLATION_$x = $_\n";
}
