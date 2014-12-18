use strict;
use Data::Dumper;

my %chr;
for (1..22){
   $chr{'chr'.$_} = '';
}
$chr{'chrX'} = '';
$chr{'chrY'} = '';
$chr{'chrM'} = '';
$chr{'chrMT'} = '';

my $encode = shift;
my $needS = shift;
#my $biomart = shift;

open EN, "$encode";
while ( <EN> ){
  chomp;
  next if /^Ensembl/;
  my @cols = split /\t/;
  my $chr = $cols[2];
  $chr = 'chr'.$chr;
  next unless (exists($chr{$chr}));
  my $start = $cols[3];
  my $end = $cols[4];
  my $strand = ($cols[5] == 1)? '+':'-';
  my $pstart = $start - $needS;
  my $pend = $start + $needS;
  if ($strand eq '-'){
    $pstart = $end - $needS;
    $pend = $end + $needS;
  }
  print "$chr\t$pstart\t$pend\t$start\t$end\t$strand\t$cols[6]\t$cols[1]\n";
}
close EN;
