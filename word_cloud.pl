use strict;

open IN, shift;
while ( <IN> ){
   chomp;
   my @cols = split /\s+/;
   foreach my $tag (@cols){
      $tag =~ /^(.+?)(\d+)$/;
      my $word = $1;
      my $count = $2;
      for (1..$count){
         print "$word ";
      }
   }
}

print "\n";
close IN;
