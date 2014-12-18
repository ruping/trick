#!/usr/bin/perl

use strict;
use Getopt::Long;

my $pattern;

GetOptions (
             "pattern|p=s" => \$pattern,
           );

my $title;
my %seq;

my $rc_pattern = reverse($pattern);
$rc_pattern =~ tr/ACGT/TGCA/;

while ( <> ){
   chomp;
   if (/^>/){
     $title = $_;
   }
   else{
     s/\n//g; 
     s/\s//g;
     $seq{$title} .= $_;
   }
}

foreach my $title (keys %seq){
   if ($seq{$title} =~ /$pattern/ or $seq{$title} =~ /$rc_pattern/){
       print "$title\n$seq{$title}\n";
   }
}


exit;
