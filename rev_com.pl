use strict;

my $seq = shift;

$seq =~ tr/ACGT/TGCA/;

$seq = reverse($seq);

print "$seq\n";
