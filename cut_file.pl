#!/usr/bin/perl -w
#cut first n line of a file

use strict;
use Getopt::Long;

my $line_s;
my $line_e;

GetOptions("ls=i" => \$line_s,
           "le=i" => \$line_e,
          );

print STDERR "$line_s\t$line_e\n";

my $line_num = 1;
while ( <> ){
  if ($line_num >= $line_s and $line_num <= $line_e) {print $_;}
  elsif ($line_num > $line_e){last;}
  $line_num++;
}
exit;
