use strict;

my $biomart = shift;
my $ensembl = shift;

open IN, "$biomart";
my %mapping;
while ( <IN> ){
  next if /^Ensembl/;
  chomp;
  my ($ens, $ensName, $ensType, $ensDes, $entrez, $wikiName, $wikiDes) = split /\t/;
  if ($entrez != ''){
    $mapping{$ens}{$entrez} = '';
  }
}
close IN;

open IN2, "$ensembl";
while ( <IN2> ){
   chomp;
   my @cols = split /\t/;
   my $ens = $cols[3];
   $ens =~ s/\.\d+$//;
   my $count = round($cols[7]);
   if (exists ($mapping{$ens})){
     foreach my $entrez (keys %{$mapping{$ens}}){
         print "$ens\t$entrez\t$count\n";
     }
   }
}
close IN2;

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}

