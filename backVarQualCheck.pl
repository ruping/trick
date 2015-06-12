use strict;
use Data::Dumper;

my $file = shift;
my $foundTh = shift;
my $ratioTh = shift;

if ($foundTh eq ''){
  $foundTh = 5;
}
if ($ratioTh eq ''){
  $ratioTh = 0.5;
}

open IN, "$file";
my %colnames;
my %colindex;
while ( <IN> ) {
  chomp;
  if (/^[\#]?chr\t/) { #it is header
    my @cols = split /\t/;
    for(my $i = 0; $i <= $#cols; $i++) {
      $colindex{$cols[$i]} = $i;
      $colnames{$i} = $cols[$i];
    }
    printf("%s\n", join("\t", @cols));
  } else {  #normal lines

    my @cols = split /\t/;
    my $colSomatic = $cols[$colindex{'somatic'}];
    my $colGermline = $cols[$colindex{'germline'}];
    my $colFounds = $cols[$colindex{'founds'}];
    my $colFreq = $cols[$colindex{'freq'}];
    my $colId = $cols[$colindex{'id'}];
    my $colFunction = $cols[$colindex{'function'}];

    if ( $colSomatic ne 'NA' ) { #check somatic ones
      my $skeep = 0;
      my @somaticSamps = split(',', $cols[$colindex{'somatic'}]);
      foreach my $somaticSamp (@somaticSamps) {
        $somaticSamp =~ /^(.+?)\[/;
        my $sSamp = $1;
        my $oriCall = $cols[$colindex{$sSamp}];
        if ($oriCall =~ /^(.+?)\|(.+?)$/) {
          my $oriMaf = $1;
          my $oriQual = $2;
          if ( $oriQual > 40 ) {
            $skeep = 1;
          } #good qual
        } #called samtools
      } #foreach somatic sample
                             #### common ####     #################### possible common #################     ### only somatic ###
      if ( $skeep == 0 and ( $colFreq > 0.005 or ($colId ne '.' and $colFreq eq 'NA' and $colFounds >= 3) or $colGermline eq 'NA' ) ) {
        next;
        #printf("%s\n", join("\t", @cols));
                                  ### common ###      #################### possible common #################
      } elsif ( $skeep == 1 and ( $colFreq > 0.01 or ($colFreq eq 'NA' and $colId ne '.' and $colFounds >= 3) ) ) {
        next;
      }
    } #check somatic ones

    if ( $colGermline ne 'NA' ) { #check germline ones
      my $gfounds = 0;
      my @germSamps = split(',', $cols[$colindex{'germline'}]);
      foreach my $gSamp (@germSamps) {
        my $oriCall = $cols[$colindex{$gSamp}];
        if ($oriCall =~ /^(.+?)\|(.+?)$/) {
          my $oriMaf = $1;
          my $oriQual = $2;
          if ( $oriQual > 30 ) {
            $gfounds += 1;
          }                     #good qual
        }                       #original call
      }                         #check germline ones
      if ( $colFounds >= $foundTh ) {
        my $gfoundRatio = sprintf("%.2f", $gfounds/$colFounds);
        if ($gfoundRatio < $ratioTh) {
          next;
          #printf("%s\n", join("\t", @cols));
        }
      }
    }  #check germline ones

    printf("%s\n", join("\t", @cols));

  } #normal lines

}
close IN;
