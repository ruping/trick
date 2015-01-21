use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Data::Dumper;

my $file = shift;

my %colnames;

open IN, "$file";
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  if ($_ =~ /^[^#]?[cC]hr\t/){
    #now it is the header
    for(my $i = 0; $i <= $#cols; $i++){
       $colnames{$cols[$i]} = $i;
    }
    print "$_\tfreq\n";
    #print STDERR Dumper(\%colnames);
  } #the clonames has been set

  else {
    #now it is the real stuff
    my $freq = -1;
    if ($cols[$colnames{'clinical'}] =~ /\;CAF\=\[([\d\.\,]+)\]\;/){
      my @freqs = split (/,/, $1);
      shift @freqs;
      #print STDERR Dumper(\@freqs);
      $freq = max(@freqs);
    } elsif ($cols[$colnames{'id'}] =~ /^[\dKGESP]+\=([\d\.]+)$/){
      $freq = $1;
    }
    print "$_\t$freq\n";
  }
}
close IN;
