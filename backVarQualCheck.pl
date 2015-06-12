use strict;
use Data::Dumper;

my $file = shift;
my $foundTh = shift;

if ($foundTh eq ''){
  $foundTh = 5;
}

open IN, "$file";
my %colnames;
my %colindex;
while ( <IN> ) {
  chomp;
  if (/^[\#]?chr\t/) {
    #it is header
    my @cols = split /\t/;
    for(my $i = 0; $i <= $#cols; $i++) {
      $colindex{$cols[$i]} = $i;
      $colnames{$i} = $cols[$i];
    }
    printf("%s\n", join("\t", @cols));
  } else {
    my @cols = split /\t/;
    if ($cols[$colindex{'somatic'}] ne 'NA') { #check somatic ones
      my $keep = 0;
      my @somaticSamps = split(',', $cols[$colindex{'somatic'}]);
      foreach my $somaticSamp (@somaticSamps) {
        $somaticSamp =~ /^(.+?)\[/;
        my $sSamp = $1;
        my $oriCall = $cols[$colindex{$sSamp}];
        if ($oriCall =~ /^(.+?)\|(.+?)$/) {
          my $oriMaf = $1;
          my $oriQual = $2;
          if ( $oriQual > 40 ) {
            $keep = 1;
          } #good qual
        } #called samtools
      } #foreach somatic sample

      if ($keep == 0){
        next;
      }
    } #somatic ones
    my @recSamps = split(',', $cols[$colindex{'germline'}]);
  }
}
close IN;
