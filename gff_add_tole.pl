use strict;

my $tole = shift;

open IN, shift;
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  $cols[3] -= $tole if ($cols[3] > $tole);
  $cols[3]  = 1 if ($cols[3] <= $tole);
  $cols[4] += $tole;
  print join("\t", @cols);
  print "\n";
}
close IN;
