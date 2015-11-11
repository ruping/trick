use strict;

my $biomart = shift;
my $ensembl = shift;
my $rpkm = shift;

open IN, "$biomart";
my %mapping;
while ( <IN> ){
  next if /^Ensembl/;
  chomp;
  my ($ens, $ensName, $ensType, $ensDes, $entrez, $wikiName, $wikiDes) = split /\t/;
  if ($entrez != ''){
    $mapping{$ens}{$entrez} = $ensName;
  }
}
close IN;

open IN2, "$ensembl";
while ( <IN2> ){
   next if /^#/;
   chomp;
   my @cols = split /\t/;
   my $ens = ($rpkm eq '')? $cols[3]:$cols[0];
   my $length = ($rpkm eq '')? $cols[8]:0;
   $ens =~ s/\.\d+$//;
   my $count = ($rpkm eq '')? round($cols[7]): $cols[1];
   if (exists ($mapping{$ens})){
     foreach my $entrez (keys %{$mapping{$ens}}){
         print "$ens\t$entrez\t$mapping{$ens}{$entrez}\t$count\t$length\n";
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

