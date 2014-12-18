#!/usr/bin/perl -w
#this is a script to get the sequence of the reference genome by entering a region

use strict;
use Getopt::Long;

my %opt = (
            'chr' => undef,
	    'start' => undef,
	    'end' => undef
          );

my $genome;

GetOptions(
            "chr=s" => \$opt{chr},
	    "start=i" => \$opt{start},
	    "end=i" => \$opt{end},
            "genome=s"=> \$genome,
	  );


my %genome = ();

open IN, "$genome";

undef $/;
$/ = '>';

while ( <IN> )
 {
   next if (/^>$/);
   s/\n>//;
   /^chr(\w+).*?\n(.*)/s;
   my $chr = $1;
   my $seq = $2;
   $seq =~ s/\n//g;
   $genome{$chr} = $seq;
 }
close IN;    #read Genome into memory

my $seq = substr ($genome{$opt{chr}}, $opt{start}-1, $opt{end}-$opt{start}+1);

print "$seq\n";

exit;
