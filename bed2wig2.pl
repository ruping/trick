use strict;
use Data::Dumper;

my $bed = shift;

my @chrs;
for (1..22){
  push(@chrs, $_);
}
push(@chrs, 'X');
push(@chrs, 'Y');
print STDERR Dumper(\@chrs);


open IN, "$bed";
my %bed;
while ( <IN> ){
  chomp;
  my ($chr, $start, $end, $tags) = split /\t/;
  push (@{$bed{$chr}{$start}}, $tags);
}
close IN;

print STDERR "all bed line stored\n";


foreach my $chr (@chrs) {
  foreach my $start (sort keys %{$bed{$chr}}) {
    print "fixedStep chrom=$chr start=$start step=100000 span=100000\n";
    foreach my $starts ( @{$bed{$chr}{$start}} ) {
      print "$starts\n";
    }
  }
}
