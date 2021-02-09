use strict;

my $hyperres = shift;
my $randomres = shift;

my %doublesig;

open IN, "$hyperres";
while ( <IN> ){
  chomp;
  my ($pathway, $id, $source, $genes, $matchgenes, $ngs, $nquery, $nmatch, $pvalue, $fdr) = split /\t/;
  next if ($ngs <= 5);
  my $name = join("\t", ($pathway, $id, $source, $genes, $matchgenes, $ngs, $nquery, $nmatch));
  $doublesig{$name}{'hyper'} = $fdr;
}
close IN;

open IN, "$randomres";
while ( <IN> ){
  chomp;
  my ($pathway, $id, $source, $genes, $matchgenes, $ngs, $nquery, $nmatch, $nbetter, $pvalue, $fdr) = split /\t/;
  my $name = join("\t", ($pathway, $id, $source, $genes, $matchgenes, $ngs, $nquery, $nmatch));
  $doublesig{$name}{'random'} = $fdr;
}
close IN;

foreach my $pname (sort {$doublesig{$a}{'random'} <=> $doublesig{$b}{'random'}} keys %doublesig) {
  my $fdrHyper = $doublesig{$pname}{'hyper'};
  next if ($doublesig{$pname}{'random'} eq '');
  my $fdrRandom = $doublesig{$pname}{'random'};
  if ($fdrHyper < 0.02 and $fdrRandom < 0.05){
    print "$pname\t$fdrHyper\t$fdrRandom\n";
  }
}
