use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Data::Dumper;

my $file = shift;
my $type = shift;  #snv or indel

if ($type eq ''){$type = 'snv'};

print STDERR "type is $type\n";

my %colnames;

open IN, "$file";
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  if ($_ =~ /^[\#]?[cC]hr\t/){
    #now it is the header
    for(my $i = 0; $i <= $#cols; $i++){
       $colnames{$cols[$i]} = $i;
    }
    print "$_\tfreq\n";
    #print STDERR Dumper(\%colnames);
  } #the clonames has been set

  else {
    #now it is the real stuff
    my $ref = $cols[$colnames{'ref'}];
    my $alt = $cols[$colnames{'alt'}];

    my $freq = -1;
    if ($cols[$colnames{'clinical'}] =~ /\;CAF\=\[([\d\.\,]+)\]\;/){
      my @freqs = split (/\,/, $1);

      $cols[$colnames{'clinical'}] =~ /\;REFALT\=([ACGT\-\,]+)$/;
      my @alleles = split(/\,/, $1);
      my $index = -1;

      for (my $i = 0; $i <= $#alleles; $i++) {
        if ($type eq 'snv'){
          if ($alleles[$i] eq $alt) {
            $index = $i;
            last;
          }
        }
      } #find index

      if ($type =~ /indel/) {
        my $knownref;
        my $knownalt;
        if ($alleles[1] =~ /^$alleles[0]/){ #insertion
          $knownref = '-';
          ($knownalt = $alleles[1]) =~ s/^$alleles[0]//;
        } elsif ($alleles[1] =~ /$alleles[0]$/) {
          $knownref = '-';
          ($knownalt = $alleles[1]) =~ s/$alleles[0]$//;
        } elsif ($alleles[0] =~ /^$alleles[1]/) {
          $knownalt = '-';
          ($knownref = $alleles[0]) =~ s/^$alleles[1]//;
        } elsif ($alleles[0] =~ /$alleles[1]$/) {
          $knownalt = '-';
          ($knownref = $alleles[0]) =~ s/$alleles[1]$//;
        }
        if ($knownref eq $ref and $knownalt eq $alt){
          $index = 1;
        } else {
          if ($type eq 'indelclean'){
            next;
          }
        }
      }

      if ($index != -1){
        $freq = $freqs[$index];
        if ($freq eq ''){
          shift @freqs;
          $freq = max(@freqs);
        }
      } else {
        shift @freqs;
        $freq = max(@freqs);
      }
    } elsif ($cols[$colnames{'id'}] =~ /^[\dKGESP]+\=([\d\.]+)$/){
      $freq = $1;
    }
    if ($type eq 'indelclean') {
      if ($freq eq $cols[$colnames{'freq'}]){
        print "$_\n";
      } else {
        next;
      }
      next;
    }
    print "$_\t$freq\n";
  }
}
close IN;
