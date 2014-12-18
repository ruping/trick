use strict;

#@SRR064286.1 HWI-EAS418:1:4:0:564 length=100
#CACCGCAAGAACCTTAGCACANAGAAAGGNNTCAACAAANNTTNGGATCCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
#+SRR064286.1 HWI-EAS418:1:4:0:564 length=100

open IN, shift;

while ( <IN> ){
  chomp;
  if ($_ =~ /^\@(\S+)\s+(\S+)\s+length\=100$/){
     my $name = $1;
     my $name1 = '@'.$name.'/1';
     my $name2 = '@'.$name.'/2';
     my $seq = <IN>;
     chomp($seq);
     my $seq1 = substr($seq, 0, 50);
     my $seq2 = substr($seq, -50);
     my $third = <IN>;
     my $third1 = '+';
     my $third2 = '+';
     my $qual = <IN>;
     chomp($qual);
     my $qual1 = substr($qual, 0, 50);
     my $qual2 = substr($qual, -50);
     print "$name1\n$seq1\n$third1\n$qual1\n";
     print STDERR "$name2\n$seq2\n$third2\n$qual2\n";
  }
}
close IN;

