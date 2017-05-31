use strict;
use File::Glob ':glob';
use Data::Dumper;

my $house = shift;
my $dir = shift;
my $out = shift;

open IN, "$house";
my @house;
while ( <IN> ) {
  chomp;
  push(@house, $_);
}
close IN;
print Dumper(\@house);

foreach my $hkg (@house) {
 my $cmd = "grep \"$hkg\$\" $dir/\*/03_STATS/\*gencode\*rpkm >>$out";
 system($cmd);
}

exit 0;

