#!/usr/bin/perl

#TODO: 1) use another file to mask the record in this file(both vcf file)
#      by the comparison of coordinates and genotypes in 2 vcf files
#      2) dbSNP id mask for vcf calls
#      written by Ruping Sun
#      Department of Systems Biology, Columbia University
#      rs3412@columbia.edu


use strict;
use Data::Dumper;

my $maskfile;
my $original;
my $zip = 0;
my $snp = 0;
my $t=0; #tolerant

my %printer;

while ($ARGV[0]) {
  my $arg = shift @ARGV;
  if ($arg =~ /-m/){$maskfile    = shift @ARGV;}
  elsif ($arg =~ /-t/){$t        = shift @ARGV;}
  elsif ($arg =~ /-o/){$original = shift @ARGV;}
  elsif ($arg =~ /-z/){$zip = 1;}
  elsif ($arg =~ /-s/){$snp = 1;}
  elsif ($arg =~ /-h/){print "useage: compare_2vcf.pl -o: origninal_filename -m: maskfile [-s:snpdb]\n";}
}

#the original file is loaded in a hash
my $oopen;
if ($zip == 1) {
  $oopen = "bgzip -dc $original |";
} else {
  $oopen = "$original";
}
open ORI, "$oopen";
my %original;
my $vcfheader = '';
while ( <ORI> ) {
  if ( /^#/ ){
    if ($snp == 0){
       next;
    } else {
       $vcfheader .= $_;
       next;
    }
  }
  chomp;
  my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $GT, $GT2) = split /\t/;
  my $chr = $CHROM;
  $chr =~ s/chr//;
  $chr = 'MT' if ($chr eq 'M');
  my $start = $POS;
  my $end = $start;
  push (@{$original{$chr}{$start}{$end}}, $_);
}
close ORI;


#the mask file is loaded in a hash as well
my $mopen;
if ($zip == 1){
  $mopen = "bgzip -dc $maskfile |";
} else {
  $mopen = "$maskfile";
}
open MASK, "$mopen";
my %mask;
while ( <MASK> ) {
  next if /^#/;
  chomp;
  my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $GT, $GT2) = split /\t/;
  my $chr = $CHROM;
  $chr =~ s/chr//;
  $chr = 'MT' if ($chr eq 'M');
  my $start = $POS;
  my $end = $start;
  push (@{$mask{$chr}{$start}{$end}}, $_);
}
close MASK;


my $counter;
foreach my $chr (sort {$a cmp $b} keys %original) {
  my @sortedpos = sort {$a <=> $b} keys %{$mask{$chr}};  #each chromosome sorted mask record into a array
  my $ptr = 0;
  foreach my $start (sort {$a <=> $b} keys %{$original{$chr}}){
    foreach my $end (sort {$a <=> $b} keys %{$original{$chr}{$start}}) {

       while (($ptr<$#sortedpos) and ((maxmaskend($sortedpos[$ptr],$chr)+$t) < $start)){$ptr++;}  #set pointer close to the mask_start
       #add 1 by 1 to test
       my $off=0;
       while (($ptr+$off<$#sortedpos) and (($sortedpos[$ptr+$off]-$t) <= $end)) {
         my $mstart = $sortedpos[$ptr+$off];
	 foreach my $mend (sort {$a <=> $b} keys %{$mask{$chr}{$mstart}}) {
	  #test if overlapped
	  if ( $end >= ($mstart-$t) and $start <= ($mend+$t) ) {
	     foreach my $mask (@{$mask{$chr}{$mstart}{$mend}}) {
	       push (@{$printer{$chr}{$start}{$end}}, $mask);
	     }
	   }
	 }
         $off++;
       }

       if ($snp == 0) {
         foreach my $record (@{$original{$chr}{$start}{$end}}) {
           print "###\t$record\t";
           if (defined $printer{$chr}{$start}{$end}) {
             print "\tfound\n";
           } else {
             print "\tmissing\n";
           }
           foreach my $printer (@{$printer{$chr}{$start}{$end}}) {
             print "$printer\n";
           }
         }
       } elsif ($snp == 1) {
         print "$vcfheader";
         foreach my $record (@{$original{$chr}{$start}{$end}}) {
           if (defined $printer{$chr}{$start}{$end}) {
             my @replaceinfo = split(/\t/, ${$printer{$chr}{$start}{$end}}[0]);
             my @record = split(/\t/, $record);
             $record[2] = $replaceinfo[2];
             $record = join("\t",@record);
             print "$record\n";
           }
         }
       }

       %printer = undef;
       $counter++;
       print STDERR "$counter\n" if ($counter%100000 == 0);

    }
  }
}


=pod
#now print the result
foreach my $chr (sort {$a cmp $b} keys %original){
  foreach my $start (sort {$a <=> $b} keys %{$original{$chr}}){
    foreach my $end (sort {$a <=> $b} keys %{$original{$chr}{$start}}){
      foreach my $record (@{$original{$chr}{$start}{$end}}) {
        print "$record\n";
	foreach my $printer (@{$printer{$chr}{$start}{$end}}) {
	  print "$printer\n";
	}
      }
    }
  }
}
=cut


sub maxmaskend{
  my ($mstart, $chr) = @_;
  my @tmp = sort {$b <=> $a} keys %{$mask{$chr}{$mstart}};
  return $tmp[0];
}
