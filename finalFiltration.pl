use strict;
use Data::Dumper;

my @name;
my %colnames;
my %found;

#my $rechecked = shift;
my $file =  shift;

#open RE, "$rechecked";
#my %re;
#while ( <RE> ) {
#   chomp;
#   my @cols = split /\t/;
#   my $coor = $cols[0].':'.$cols[1];
#   $re{$coor}{'mean'} = $cols[$#cols-1];
#   $re{$coor}{'median'} = $cols[$#cols];
#}
#close RE;

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
    my $segdupScore = 0;
    #if (exists ($re{$coor})) {
    #  $cmean = $re{$coor}{'mean'};
    #  $cmedian = $re{$coor}{'median'};
    #} else {
    #  print STDERR "not found consecutive mismatch mean and median for $coor\n";
    #}
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
    }                           #iterator
    if ($junc < 4.1) {
      print "$_\t$junc\n";
    }
  }                             #else
}                               #while
close IN;
