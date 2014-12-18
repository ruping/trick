use strict;

my $genome = shift;

open HS, "$genome";
my $chr = undef;
my $countAll = 0;
my $countN = 0;
while ( <HS> ) {
  if (/^>(\w+).*?\n$/) {
     $chr = $1;
     $chr =~ s/^chr//;
  }
  else {
     s/\n//g;
     s/\s//g;
     my $length = length($_);
     $countAll += $length;
     if ($_ =~ /NNNNNNNN/){
        $countN += $length;
     }
  }
}
close HS;
my $rest = $countAll-$countN;
print "$countAll\t$countN\t$rest\n";
