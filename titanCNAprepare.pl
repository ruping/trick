use strict;
use Data::Dumper;
use File::Glob ':glob';
use File::Basename;


my $data = shift;
my $somaticInfo = shift;
my $pairedCall = shift;
my $split = 1;

my %somatic;
my %germline;  #may have multiple tumors
if ($somaticInfo ne '' and -s "$somaticInfo") {

  open IN, "$somaticInfo";
  while ( <IN> ){
    chomp;
    s/[\s\n]$//;
    my @columns = split /\t/;
    my $tumor = $columns[0];
    my $normal = $columns[1];

    $somatic{$tumor} = $normal;
    push(@{$germline{$normal}}, $tumor) if $normal ne 'undef';
  }
  close IN;
  print STDERR Dumper (\%somatic);
  print STDERR Dumper (\%germline);
}


my $outdir = "./titan/";
unless (-e "$outdir"){
  system("mkdir -p $outdir");
}


if ($split == 1) {

  my %colnames;
  my %colindex;
  my %fhs;
  open IN, "$data";
  while (<IN>) {
    chomp;
    my @cols = split /\t/;
    if ($_ =~ /^[\#]?chr\t/) {
      $_ =~ s/^\#//;
      for (my $i = 0; $i <= $#cols; $i++) {
        $colnames{$i} = $cols[$i];
        $colindex{$cols[$i]} = $i;
      }
      next;
    } else {
      my $chr = $cols[$colindex{'chr'}];
      my $pos = $cols[$colindex{'pos'}];
      my $ref = $cols[$colindex{'ref'}];
      my $alt = $cols[$colindex{'alt'}];
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colnames{$i} =~ /^(.+?)maf$/){  #now it is sample maf

          my $sample = $1;

          if (exists($germline{$sample})) {  #it is a blood
            my $calledBlood = $cols[$i-1];
            if ( $pairedCall == 1 ) {
              $calledBlood = $cols[$colindex{${$germline{$sample}}[0]}];
            }
            if ($calledBlood =~ /\|/) {   #originally called
              if ($cols[$i] =~ /\|/) { #split the var surrounding information
                my @infos = split(/\|/, $cols[$i]);
                my $bmaf = $infos[0];
                my $bendsratio = $infos[1];
                my ($bcmean, $bcmedian) = split(',', $infos[2]);
                my $strandRatio = $infos[3];
                my $badQualFrac = $infos[4];
                if ($bendsratio <= 0.9 and ($strandRatio != 0 and $strandRatio != 1) and $badQualFrac < 0.6 and (($bcmean+$bcmedian) < 5.5 or $bcmedian <= 2)) { #likely true event

                  foreach my $tumorSamp (@{$germline{$sample}}) {   ##now should start checking for each tumor samples

                    my $indexts = $colindex{$tumorSamp.'maf'};
                    if ($cols[$indexts] =~ /\|/) { #split the var surrounding information
                      my @tsinfo = split(/\|/, $cols[$indexts]);
                      my $tsmaf = $tsinfo[0];
                      my $tsendsratio = $tsinfo[1];
                      my ($tscmean, $tscmedian) = split(',', $tsinfo[2]);
                      my $tsd = $cols[$indexts+1];
                      if ($tsendsratio <= 0.9 and (($tscmean+$tscmedian) < 5.5 or $tscmedian <= 2)) {  #likely true event, start printing
                        my $fh = $tumorSamp;
                        unless (-e "$outdir/$tumorSamp\_titan") {
                          open ( my $fh, ">>", "$outdir/$tumorSamp\_titan" )  || die $!;
                          $fhs{$tumorSamp} = $fh;
                          print {$fhs{$tumorSamp}} "chr\tpos\tref\trefCount\talt\taltCount\n";
                        }
                        my $NrefCount = 0;
                        my $refCount = 0;
                        $NrefCount = round($tsmaf*$tsd);
                        $refCount = $tsd - $NrefCount;
                        if (($refCount + $NrefCount) >= 5) {
                          print {$fhs{$tumorSamp}} "$chr\t$pos\t$ref\t$refCount\t$alt\t$NrefCount\n";
                        }
                      } #true event print
                    } #split tumor info
                  }  ##now should start checking for each tumor samples

                } #true blood event
              } #split blood recheck info
            } #originally called
          } #blood

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
