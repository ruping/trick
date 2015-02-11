use strict;

#my $num = shift;
#my $int = round($num);
#print "$int\n";
my $blood = shift;
my $tumor = shift;
my 


open IN, "$blood";
while ( <IN> ){
  chomp;
  if (/^#/){
    print "$_\n";
    next;
  }
  my ($chr, $pos, $id, $ref, $alt, $bloodmaf, $bloodd, $cmeanav, $cmedianav) = split /\t/;
  my $vard = &round($bloodmaf*$bloodd);
  next if ($cmeanav >= 3 and $cmedianav >= 3);
  if ($bloodd >= 7 and $bloodmaf >= 0.25 and $bloodmaf <= 0.75){
    print "$_\n";
  }
}
close IN;




sub round {
  my $number = shift;
  my $tmp = int($number);
  if ($number >= ($tmp+0.5)){
    $tmp++;
  }
  return $tmp;
}
