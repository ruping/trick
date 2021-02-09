use strict;
use Data::Dumper;

my $file = shift;

my %data;

open IN, "$file";
my %colnames;
my %colindex;
my @samples;
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  if ($_ =~ /^gene/) {  #header
    for (my $i = 0; $i <= $#cols; $i++) {
      $colnames{$i} = $cols[$i];
      $colindex{$cols[$i]} = $i;
    }
  } #header
  else { #normal lines
    my $gene = $cols[$colindex{'gene'}];
    #print "$gene\n";
    for (my $i = 1; $i <= $#cols; $i++) {
      my $sample = $colnames{$i};
      if ( scalar(@samples) < 21 ) {
        push(@samples, $sample);
      }
      $data{$gene}{$sample} = $cols[$i];
    }
  }
}
close IN;

#print STDERR Dumper(\%colnames);
#print STDERR Dumper(\%colindex);
#print Dumper(\%data);
#print Dumper(\@samples);

srand();
my $b = 1;
while ($b <= 1000) {
  my %rinds;
  my $s = 1;  #7 samples
  while ( $s <= 7 ) {
    my $rind = int(rand(21));
    while ( exists($rinds{$rind}) ) {
      $rind = int(rand(21));
    }
    $rinds{$rind} = '';
    $s++;
  } #7 samples
  #print Dumper(\%rinds);

  #now the samples are chose, do genes
  foreach my $srindp (sort keys %rinds) {  #samples
    my $samplep = $samples[$srindp];
    print "b$b\t$samplep\n";
  }                             #each sample

  open OUT, ">b$b";
  foreach my $gene (keys %data) {   #each gene
    my $found  = 0;
    foreach my $srind (sort keys %rinds) {  #samples
      my $sample = $samples[$srind];
      if ($data{$gene}{$sample} != 0){
        $found = 1;
      }
    } #each sample
    if ($found != 0){
      print OUT "$gene\n";
    }
  }
  close OUT;

  #now run pathway
  my $cmd = "perl /tools/trick/runPathwayEnrichment.pl b$b /cygdrive/h/annotation/CPDB_pathways_genes.tab";
  RunCommand($cmd, 0 , 0);

  print STDERR "b$b\n";

  $b++;
}

exit 0;


sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}
