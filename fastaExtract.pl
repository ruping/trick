use strict;
use Data::Dumper;

my $genome = shift;
my $chrs = shift;

my @chrs = split(',', $chrs);
my %loci;
foreach my $loci (@chrs) {
  $loci =~ /^(chr[0-9XYMT]+)\:(\d+)\-(\d+)$/;
  my $chr = $1;
  my $start = $2;
  my $end = $3;
  $loci{$chr}{'start'} = $start;
  $loci{$chr}{'end'} = $end;
}
print STDERR Dumper(\%loci);

open HS, "$genome";
my %genome = ();
my $chr = undef;
while ( <HS> ) {
  if (/^>(\w+).*?\n$/) {
    $chr = $1;
  }
  else {
    s/\n//g;
    s/\s//g;
    $genome{$chr}.=$_;
  }
}
close HS;
print STDERR "genome loaded\n";


foreach my $chr (keys %loci) {

  my $start = $loci{$chr}{'start'};
  my $end = $loci{$chr}{'end'};
  my $seq = substr($genome{$chr}, $start-1, $end-$start+1);

  my $identifier = $chr.'_'.$start.'_'.$end;
  print ">$identifier\n";
  print "$seq\n";

}

exit 0;
