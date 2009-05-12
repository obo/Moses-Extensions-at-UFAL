#!/usr/bin/perl

BEGIN { use Cwd qw/ abs_path /; use File::Basename; $script_dir = dirname(abs_path($0)); push @INC, "$script_dir/../perllib"; }
use RegTestUtils;

$x=0;
while (<>) {
  chomp;

  if (/^Finished loading LanguageModels/) {
    my $time = RegTestUtils::readTime($_);
    print "LMLOAD_TIME ~ $time\n";
  }
  if (/^Finished loading phrase tables/) {
    my $time = RegTestUtils::readTime($_);
    print "PTLOAD_TIME ~ $time\n";
  }
  next unless /^BEST TRANSLATION:/;
  my $pscore = RegTestUtils::readHypoScore($_);
  print "SCORE_$x = $pscore\n";
  $x++;
}
