use strict;
use Data::Dumper;

my $mindy = shift;
my $mut = shift;
my $reverse = shift;

open MIN, "gzip -dc $mindy |";
my %mindy;
while ( <MIN> ){
  chomp;
  next if ($_ =~ /^modEntrez\t/);
  my ($modEntrez, $modSymbol, $tfEntrez, $tfSymbol, $nrtargets, $p, $CINDYTargetGenes) = split /\t/;
  if ($reverse eq ''){
    $mindy{$tfSymbol}{$modSymbol} = 0;
  } else {
    $mindy{$modSymbol}{$tfSymbol} = 0;
  }
}
close MIN;

#print Dumper(\%mindy);

open MUT, "$mut";
my %mut;
while ( <MUT> ) {
  next if ($_ =~ /^[#@]/ or $_ =~ /^chr\t/);
  my @cols = split /\t/;
  my $gene = $cols[0];
  $mut{$gene} = '';
}
close MUT;


if ($reverse eq ''){
  foreach my $tfSymbol (keys %mindy){
    my $modN = 0;
    my $mutN = 0;
    foreach my $modSymbol (keys %{$mindy{$tfSymbol}}){
      $modN ++;
      if (exists $mut{$modSymbol}){
         $mindy{$tfSymbol}{$modSymbol} ++;
         $mutN ++;
      }
    }
    my $ratio = sprintf("%.3f", $mutN/$modN);
    print "$tfSymbol\t$modN\t$mutN\t$ratio\n";
  }
} else {
  foreach my $modSymbol (keys %mindy){
    my $tfN = 0;
    my $mutN = 0;
    if (exists $mut{$modSymbol}) {
       $tfN = 1;
    }
    foreach my $tfSymbol (keys %{$mindy{$modSymbol}}){
      $mutN ++;
    }
    print "$modSymbol\t$tfN\t$mutN\n";
  }
}
