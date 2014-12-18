use strict;

#bin    score   tName   tSize   tStart  tEnd    qName   qSize   qStrand qStart  qEnd    id
#607     32612   chr1    197195432       3000001 3001179 chr3    159599783       -       134441891       134442892       5003118

open IN, shift;

my $line = 1;

while ( <IN> ) {
   next if /^#/;
   chomp;
   my ($bin, $score, $tName, $tSize, $tStart, $tEnd, $qName, $qSize, $qStrand, $qStart, $qEnd, $id) = split /\t/;
   my $chr1 = $tName;
   my $chr2 = $qName;
   my $source = 'ucsc';
   my $type = 'selfChain';
   my $start1 = $tStart;
   my $start2 = $qStart;
   my $end1 = $tEnd;
   my $end2 = $qEnd;
   my $score = '.';
   my $strand = '.';
   my $phase = '.';
   my $tag = "L=$line";
   $line++;
   printf("%s\n", join("\t", $chr1, $source, $type, $start1, $end1, $score, $strand, $phase, $tag));
   printf("%s\n", join("\t", $chr2, $source, $type, $start2, $end2, $score, $strand, $phase, $tag));
}
close IN;
