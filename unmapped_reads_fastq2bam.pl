#!/usr/bin/perl
use warnings;
use strict;

my ($fastq,$sam,$outfile) = @ARGV;

unless ($outfile) {
  die "Usage is filter_unmapped_reads.pl [FastQ file] [SAM File] [File for unmapped reads]\n";
}

if (-e $outfile) {
  die "Won't overwrite an existing file, delete it first!";
}

open (FASTQ,$fastq) or die "Can't open fastq file: $!";
open (SAM,$sam) or die "Can't open SAM file: $!";
open (OUT,'>',$outfile) or die "Can't write to $outfile: $!";

my $ids = read_ids();

filter_fastq($ids);

close OUT or die "Can't write to $outfile: $!";


sub filter_fastq {

  warn "Filtering FastQ file\n";

  my ($ids) = @_;

  while (<FASTQ>) {

    if (/^@(.+?)[\/\s]/) {

      my $id = $1;
      my $seq = <FASTQ>;
      my $id2 = <FASTQ>;
      my $qual = <FASTQ>;


      unless (exists $ids->{$id}) {
	print OUT $_,$seq,$id2,$qual;
      }
    }
    else {
      warn "Line '$_' should have been an id line, but wasn't\n";
    }

  }

}


sub read_ids {

  warn "Collecting mapped ids\n";

  my $ids;

  while (<SAM>) {

    next if (/^@/);
    my ($id) = split(/\t/);
    $ids->{$id} = 1;
  }

  return $ids;
}
