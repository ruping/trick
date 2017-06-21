use strict;

my $genome = shift;
my $file = shift;

open HS, "$genome";
my %genome = ();
my $chr = undef;
while ( <HS> ) {
  if (/^>(\w+).*?\n$/) {
    $chr = $1;
    #$chr =~ s/^chr//;
  } else {
    s/\n//g;
    s/\s//g;
    $genome{$chr} += length($_);
  }
}
close HS;
print STDERR "genome loaded\n";

foreach my $chr (sort keys %genome) {
  print "$chr\t$genome{$chr}\n";
}

close IN;
