use strict;


#rep mask
my %repeatmask;
my @rs;      #start array for each chr
my $old_chr; #checking the chr
my $ptr;     #pointer for repeatmask sub
open REP, "/scratch/ngsvin/TargetSequencing/TESTING/ruping/workspace/RepeatMask/Hs.UCSC_repeats.gff";
  while (<REP>){
    next if /^#/;
    chomp;
    my @tmp = split /\t/;
    #(my $chr = $tmp[0]) =~ s/^chr(\w+)$/$1/;
    my $chr = $tmp[0];
    my $repeat_s = $tmp[3];
    my $repeat_e = $tmp[4];
    $repeatmask{$chr}{$repeat_s} = $repeat_e;
  }
close REP;
print STDERR "rep loaded\n";


open IN, shift;
while ( <IN> ){
   chomp;
   my @cols = split /\t/;
   if (repeatmask($cols[0], $cols[3]) == 0){
      print "$_".";rep=no\n";
   }
   else {
      print "$_".";rep=yes\n";
   }
}
close IN;



sub repeatmask {
    my ($chr, $coor) = @_;
    my $flag = 0;
    if ($chr ne $old_chr){
      @rs = sort {$a <=> $b} keys %{$repeatmask{$chr}};
      $ptr = 0;
    }
    while (($ptr<=$#rs) and ($repeatmask{$chr}{$rs[$ptr]} < $coor)){
      $ptr++;
    }
    if ($rs[$ptr] <= $coor){
      $flag = 1;
    }
    $old_chr = $chr;
    return $flag;
}
