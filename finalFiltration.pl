use strict;
use Data::Dumper;

my @name;
my %colnames;
my %found;

my $file = shift;
my $type = "snv";



open IN, "$file";
while ( <IN> ) {
  chomp;
  my @cols = split /\t/;
  if ($_ =~ /^[\#]?chr\t/) {
    @name = @cols;
    for (my $i = 0; $i <= $#name; $i++){
      $colnames{$name[$i]} = $i;
    }
    print "$_\tscore\n";
    next;
  } else {
    my $junc = 0;
    my $coor = $cols[0].':'.$cols[1];
    my $cmean = $cols[$colnames{'cmeanav'}];                                            #consecutive mismatch mean
    my $cmedian = $cols[$colnames{'cmedianav'}];                                        #consecutive mismatch median
    my $ref = $cols[$colnames{'ref'}];
    my $alt = $cols[$colnames{'alt'}];
    if ($ref !~ /^[ACGTN]$/ or $alt !~ /^[ACGTN]$/){                             #indel
      $type = "indel";
    }
    my $rep = $cols[$colnames{'rep'}];
    my $sc = $cols[$colnames{'sc'}];
    my $segdupScore = 0;

    unless ($cmean < 2 and $cmedian < 2) {
       $junc += ($cmean + $cmedian)/2;                                           #good ones
    }
    for (my $i = 0; $i <= $#cols; $i++) {
      if ($name[$i] =~ /func/){
        if ($cols[$i] =~ /seg(d)?up\.score\=([\d\.]+)/) {
           $segdupScore = $2;
        }
      } elsif ($name[$i] =~ /\.bam\.out\.bad/) {
          $junc += 4*$cols[$i];
      } elsif ($name[$i] eq 'rep') {
          $junc += 2*$cols[$i];
      } elsif ($name[$i] eq 'sc') {
          $junc += (($segdupScore + $cols[$i]) > 0) ? 2:0;
      }
    } #iterator
    if ($type eq 'indel' and $rep == 0 and $sc == 0 and $segdupScore == 0) {
      if  ($junc < 4.3) {
        print "$_\t$junc\n";
      }
    } else {
      if  ($junc < 4) {
        print "$_\t$junc\n";
      }
    }
  }                             #else
}                               #while
close IN;
