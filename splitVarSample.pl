use strict;
use Data::Dumper;

my @name;
my %found;

open IN, shift;
while ( <IN> ) {
  chomp;
  next if /^@/;
  my @cols = split /\t/;
  if ($_ =~ /^#chr\t/) {
    @name = @cols;
    next;
  } else {
    my %max;
    for (my $i = 0; $i <= $#cols; $i++) {
      if ($name[$i] =~ /^(AC\d+)$/){
        my $cus = $1;
        next if $cus eq 'AC19';
        next if $cus eq 'AC20';
        next if $cus eq 'AC21';
        #$max{$cus}{'called'} = $cols[$i];
        $max{$cus}{'rescue'} = $cols[$i];
        $max{$cus}{'depth'}  = $cols[$i+1];
      } elsif ($name[$i] eq 'AC3T1') {
        my $cus = "AC57";
        $max{$cus}{'called'} = $cols[$i];
        $max{$cus}{'rescue'} = $cols[$i];
        $max{$cus}{'depth'}  = $cols[$i+4];
      } elsif ($name[$i] eq 'AC3T2') {
        my $cus = "AC583";
        $max{$cus}{'called'} = $cols[$i];
        $max{$cus}{'rescue'} = $cols[$i];
        $max{$cus}{'depth'}  = $cols[$i+3];
      }
    }                           #iterator
    #find the best sample
    my @bestSamples = sort {
                            #my $ca=$max{$a}{'called'}; my $cb=$max{$b}{'called'};
                            my $ra=$max{$a}{'rescue'}; my $rb=$max{$b}{'rescue'};
                            my $da=$max{$a}{'depth'}; my $db=$max{$b}{'depth'};
                            #$cb <=> $ca ||
                            $rb <=> $ra || $db <=> $da
                           } keys %max;
    my $bestSample = $bestSamples[0];
    open OUT, ">>consecutiveMismatch/snv/$bestSample";
    print OUT "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\n";
    close OUT;
  }                             #else
}                               #while
close IN;
