use strict;

my $file = shift;

#P-values for strand bias, baseQ bias, mapQ bias and tail distance bias

open IN, "$file";
my $call;
while  ( <IN> ){
  chomp;
  if ($_ =~ /^###/){
    my @cols = split /\t/;
    $call = $cols[1]."\t".$cols[2];
  } else {
    my @cols = split /\t/;
    $cols[7] =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/;
    my $strandb = $1;
    my $baseqb = $2;
    my $mapqb = $3;
    my $tailb = $4;
    my @bias;
    push (@bias, $strandb);
    push (@bias, $baseqb);
    push (@bias, $mapqb);
    push (@bias, $tailb);
    my $bad = 0;
    if ($cols[6] < 90){
      $bad = 1;
    }
    foreach my $bias (@bias) {
      if ($bias =~ /e/) {
        $bad = 1;
      } elsif ($bias <= 0.05) {
        $bad = 1;
      }
    }
    print "$call\t$bad\n";
  }
}
close IN;
