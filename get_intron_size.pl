use strict;

open IN, shift;

while ( <IN> ) {
  chomp;
  my @cols = split /\t/;
  next if $cols[9] == 1;

  my @blockSizes = split /\,/, $cols[10];
  my @blockStarts = split /\,/, $cols[11];
  for (0..($#blockSizes-1)){
     my $intronSize = $blockStarts[$_+1]-$blockStarts[$_]-$blockSizes[$_];
     print "$intronSize\n";
  }

}

close IN;
