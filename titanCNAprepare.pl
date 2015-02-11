use strict;

#my $num = shift;
#my $int = round($num);
#print "$int\n";
my $data = shift;
my $split = shift;

my $outdir = "./";

unless ($split == 1){
  open IN, "$data";
  while ( <IN> ) {
    chomp;
    if (/^#/) {
      print "$_\n";
      next;
    }
    my ($chr, $pos, $id, $ref, $alt, $bloodmaf, $bloodd, $cmeanav, $cmedianav) = split /\t/;
    my $vard = &round($bloodmaf*$bloodd);
    next if ($cmeanav >= 3 and $cmedianav >= 3);
    if ($bloodd >= 7 and $bloodmaf >= 0.25 and $bloodmaf <= 0.75) {
      print "$_\n";
    }
  }
  close IN;
}

if ($split == 1){
  my %colnames;
  open IN, "$data";
  while (<IN>){
    chomp;
    my @cols = split /\t/;
    if ($_ =~ /^[\#]?chr\t/) {
      $_ =~ s/^\#//;
      for (my $i = 0; $i <= $#cols; $i++) {
        $colnames{$i} = $cols[$i];
      }
      next;
    } else {
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colnames{$i} =~ /^(.+?)maf$/){  #now it is maf
          my $sample = $1;
          open OUT, ">>$outdir/$sample\_titan";
          my $refCount = 0;
          my $NrefCount = 0;
          $refCount = round($cols[$i]*$cols[$i+1]);
          $NrefCount = $cols[$i+1] - $refCount;
          if (($refCount +$NrefCount) >= 3) {
            print OUT "$cols[$colnames{'chr'}]\t$cols[$colnames{'pos'}]\t$cols[$colnames{'ref'}]\t$refCount\t$cols[$colnames{'alt'}]\t$NrefCount\n";
          }
          close OUT;
        } #maf
      } #each col
    } #each non header
  } #each line
  close IN;
} #split samples


sub round {
  my $number = shift;
  my $tmp = int($number);
  if ($number >= ($tmp+0.5)){
    $tmp++;
  }
  return $tmp;
}
