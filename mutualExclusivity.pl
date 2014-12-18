use strict;

my @ileum = qw(AC444 AC445 AC446 AC516 AC517 AC518 AC519);
my %ileum;
foreach my $il (@ileum){
  $ileum{$il} = '';
}
my @rectum = qw(AC3 AC439 AC440 AC441 AC443 AC447 AC525 AC526 AC527 AC528 AC529 AC530 AC531 AC532 AC533 AC546 AC548 AC580);
my %rectum;
foreach my $rec (@rectum){
  $rectum{$rec} = '';
}

open IN, shift;
my %data;
while ( <IN> ){
   chomp;
   my ($gene, $samples, $rec, $il) = split /\t/;
   next if $gene eq '';
   my @samples = split (/\,/, $samples);
   foreach my $sample (@samples){
     $data{$gene} = \@samples;
   }
}
close IN;
print STDERR "data loaded\n";

foreach my $gene (keys %data){
  foreach my $gene2 (keys %data){
    next if $gene eq $gene2;
    my $il;
    my $rec;
    my %tmp;
    foreach my $sample (@{$data{$gene}}){
      $tmp{$sample} = '';
    }
    foreach my $sample (@{$data{$gene2}}){
      $tmp{$sample} = '';
    }
    foreach my $sample (keys %tmp){
      if (exists($ileum{$sample})){
        $il++;
      } elsif (exists($rectum{$sample})){
        $rec++;
      }
    }
    print "$gene\t$gene2\t$rec\t$il\n";
  }
}
