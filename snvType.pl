use strict;

my %mutype;

my $genome = shift;
my $mut = shift;


open HS, "$genome";
my %genome = ();
my $chr = undef;
while ( <HS> ) {
  if (/^>(\w+).*?\n$/) {$chr = $1; $chr =~ s/^chr//;}
  else { s/\n//g; s/\s//g; $genome{$chr}.=$_;}
}
close HS;
print STDERR "genome loaded\n";


open IN, "$mut";

while ( <IN> ) {
  chomp;
  next if /^[@#]/;
  next if /^chr\t/;
  my @cols = split /\t/;
  my $chr = $cols[0];
  my $pos = $cols[1];
  my $ref = $cols[3];
  my $alt = $cols[4];
  my $seq = substr($genome{$chr},$pos-2,3);
  my $realref = substr($seq, 1, 1);
  if ($realref ne $ref) {
     print STDERR "shit\t$chr\t$pos\t$ref\t$realref\n";
     exit 22;
  }
  my $oldseq = $seq;
  my $newseq = substr($seq, 1, 1, $alt);
  $newseq = $seq;

  if ($ref =~ /[GT]/){ #do reverse complementary
     $oldseq = reverse($oldseq);
     $oldseq =~ tr/ACGT/TGCA/;
     $newseq = reverse($newseq);
     $newseq =~ tr/ACGT/TGCA/;
  }

  $mutype{$oldseq."->".$newseq} ++;
}

close IN;


foreach my $mutype (sort {$mutype{$b} <=> $mutype{$a}} keys %mutype){
  print "$mutype\t$mutype{$mutype}\n";
}
