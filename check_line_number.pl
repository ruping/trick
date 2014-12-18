#!/usr/bin/perl

use strict;

my $file = shift;
my $pattern = shift;

if ($file =~ /\.gz$/){
  open IN,"gzip -d -c $file|";
} elsif ($file =~ /\.bz2$/){
  open IN,"bzip2 -d -c $file|";
} else {
  open IN,"$file";
}

my $line = 1;
while ( <IN> ){
  if ($_ =~ /$pattern/){
    print "$line\n";
  }
  $line++;
}
close IN;
$line--;
print "$line\n";
