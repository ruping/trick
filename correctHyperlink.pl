use strict;

open IN, shift;
my @colnames;
while ( <IN> ){
   chomp;
   if (/^chr/){
     @colnames = split /\t/;
     print "$_\n";
   } else {
     my @cols = split /\t/;
     for (my $i = 0; $i <= $#cols; $i++){
       if ($colnames[$i] eq 'link'){
         $cols[$i] =~ s/HYPERLINK\(/HYPERLINK\(\"/;
         $cols[$i] =~ s/\, UCSC\)/\"\, \"UCSC\"\)/;
       }
     }
     printf("%s\n", join("\t", @cols));
   }
}
close IN;
