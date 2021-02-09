use strict;

my $res = shift;

my %epi;
$epi{'E101'} = "Rectal_Mucosa";
$epi{'E102'} = "Rectal_Mucosa2";
$epi{'E103'} = "Rectal_Smooth_Muscle";
$epi{'E077'} = "Duodenum_Mucosa";
$epi{'E096'} = "Lung";
$epi{'E087'} = "Pancreatic_Islets";
$epi{'E075'} = "Colonic_Mucosa";

open IN, "$res";
my %colnames;
my %colindex;
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  for(my $i = 0; $i <= $#cols; $i++){ #clean it
    $cols[$i] =~ s/\"//g;
    $cols[$i] =~ s/[\s\n]$//;
  }
  if ($_ =~ /^#chr/) {
    for(my $i = 0; $i <= $#cols; $i++){
      if ($cols[$i] =~ /broadPeak\.sorted\.name/){
        $cols[$i] =~ /(E\d+)\-(\w+)\./;
        my $tissue = $epi{$1};
        my $mark = $2;
        $cols[$i] = $tissue.'-'.$mark;
      }
      $colindex{$cols[$i]} = $i;
      $colnames{$i} = $cols[$i];
    }
    printf("%s\n", join("\t", @cols, 'epi'));
  } else {
    my $status = 'NA';
    for($colindex{'Rectal_Mucosa-H3K27ac'}..$colindex{'Colonic_Mucosa-H3K9me3'}){
      my $idx = $_;
      my $hmark = $colnames{$idx};
      if ($cols[$idx] =~ /^Rank\_(\d+)$/){
        if ($1 < 20000){
          $status .= "$hmark\=$1," if $status ne 'NA';
          $status = "$hmark\=$1," if $status eq 'NA';
        }
      }
    }
    printf("%s\n", join("\t", @cols, $status)) if $status =~ /,\w/;
  }
}
close IN;
