use strict;

my $genome = shift;
my $file = shift;

open HS, "$genome";
my %genome = ();
my $chr = undef;
while ( <HS> ) {
  if (/^>(\w+).*?\n$/) {$chr = $1; $chr =~ s/^chr//;}
  else { s/\n//g; s/\s//g; $genome{$chr}.=$_;}
}
close HS;
print STDERR "genome loaded\n";

open IN, "$file";
while ( <IN> ){
  chomp;
  next if /^[\#]?chr\t/;
  my @cols = split /\t/;
  (my $chrom = $cols[0]) =~ s/^chr//;
  my $start = $cols[1];
  my $end   = $cols[2];
  my $id = $chrom.':'.$start.'-'.$end;
  my $seq = substr($genome{$chrom}, ($start-1), ($end-$start+1));
  print ">$id\n$seq\n";
}
close IN;
