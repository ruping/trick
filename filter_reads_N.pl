#!/usr/bin/perl

use strict;

my $title;

open IN, "$ARGV[0]";
while ( <IN> ){
   chomp;
   if (/^>/){
     $title = $_;
   }
   else{
     s/\n//g; 
     s/\s//g;
     my $numN;
     while ($_ =~ /N/g){
       $numN++;
     }
     if ($numN < 3){
       print "$title\n$_\n";
     }
   }
}
close IN;

exit;