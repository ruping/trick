use strict;

open IN, shift;
my @colnames;
while ( <IN> ){
   chomp;
   if ($_ =~ /^[#]?chr\t/){
     print "$_\n";
     $_ =~ s/^#//;
     @colnames = split /\t/;
   } else {
     my @cols = split /\t/;
     my $expressed = 0;
     my $total = 0;
     for (my $i = 0; $i <= $#cols; $i++){
       if ($colnames[$i] =~ /AC\d+R\.AC\d+/){
         $total++;
         if ($cols[$i] >= 0.1 and $cols[$i+1] >= 5){
            $expressed++;
         }
       }
     }
     if ($expressed/$total >= 0.5){
       print "$_\n";
     }
   }
}
close IN;
