use strict;

open IN, shift;
my $gene;
my $n = 0;
my %meth;
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  if ($_ =~ /^#/) {
    if ($gene ne '') {
      if ($n > 0) {   #print the old one
        my $weighted;
        my $alldepth;
        foreach my $CpG (@{$meth{$gene}}) {
           $alldepth += $CpG->{'depth'};
        }
        foreach my $CpG (@{$meth{$gene}}) {
           my $depth = $CpG->{'depth'};
           my $methy = $CpG->{'methy'};
           $weighted += ($depth/$alldepth)*$methy;
        }
        $weighted = sprintf("%.8f", $weighted);
        print "$gene\t$weighted\n";
      } else {
        print "$gene\t-1\n";
      }
    }
    ($gene = $_) =~ s/^[\#]+\t//;
    $n = 0;
  } else {
    $cols[3] =~ /^\'(\d+)\/(\d+)\'$/;
    my $methD = $1;
    my $totD = $2;
    next if $totD < 3;                  #skip exteremly low depth site
    $totD = ($totD <= 25)? $totD : 25;  #upper bound of depth
    my %tmp;
    $tmp{'depth'} = $totD;
    $tmp{'methy'} = $cols[4]/1000;
    push(@{$meth{$gene}}, \%tmp);
    $n++;
  }
}
close IN;
