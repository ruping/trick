use strict;

open IN, shift;
print "#gff3 ucsc repeatmasker track\n";

while ( <IN> ) {
   chomp;
   next if /^#/;

   my ($bin,$swScore,$milliDiv,$milliDel,$milliIns,$genoName,$genoStart,$genoEnd,$genoLeft,$strand,$repName,$repClass,$repFamily,$repStart,$repEnd,$repLeft,$id) = split /\t/;

   my $chr = $genoName;
   my $source = "ucsc";
   my $type = "repeat_region";
   my $start = $genoStart;
   my $end = $genoEnd;
   my $score = '.';
   my $tag = 'repName='.$repName.';repClass='.$repClass.';repFamily='.$repFamily;
   printf("%s\n", join("\t", $chr,$source,$type,$start,$end,$score,$strand,$score,$tag));

}

close IN;
