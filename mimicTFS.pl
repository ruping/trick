use strict;

my $sites = shift;
my $allTFs = shift;

open IN, "$sites";
my %sites;
while ( <IN> ){
  chomp;
  next if /^[\#]?chr\t/;
  my @cols = split /\t/;
  my $coor = $cols[0].':'.$cols[1].'-'.$cols[2];
  $sites{$coor} = '';
}
close IN;

open TF, "$allTFs";
my %tfs;
while ( <TF> ){
  chomp;
  s/[\s\n]//g;
  my $tf = $_;
  $tfs{$tf} = '';
}
close TF;

foreach my $coor (keys %sites){
  print "$coor";
  foreach my $tf (keys %tfs){
    print "\t$tf\=1";
  }
  print "\n";
}
