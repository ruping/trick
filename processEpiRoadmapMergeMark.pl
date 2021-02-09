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
      if ($cols[$i] =~ /broadPeak\.sorted/){
        $cols[$i] =~ /(E\d+)\-(\w+)\./;
        my $tissue = $epi{$1};
        my $mark = $2;
        $cols[$i] = $tissue.'-'.$mark;
      }
      $colindex{$cols[$i]} = $i;
      $colnames{$i} = $cols[$i];
    }
    printf("%s\n", join("\t", $cols[0], $cols[1], $cols[2], 'mark'));
  } else {
    my $status = 'NA';
    for($colindex{'Rectal_Mucosa-H3K27ac'}..$colindex{'wgEncodeRegDnaseClusteredV2.sorted'}){
      my $idx = $_;
      my $hmark = $colnames{$idx};
      if ($cols[$idx] == 1){
        $status .= "$hmark," if $status ne 'NA';
        $status = "$hmark," if $status eq 'NA';
      }
    }
    printf("%s\n", join("\t", $cols[0], $cols[1], $cols[2], $status));
  }
}
close IN;
