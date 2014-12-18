#!/usr/bin/perl

use strict;

my $seq = shift;

$seq = reverse($seq);

$seq =~ tr/ACGTacgt/TGCAtgca/;

print "$seq\n";
