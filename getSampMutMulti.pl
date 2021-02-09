#this script is for preparing multi sample mut table#
#          either for somatic or germline           #

use strict;
use Data::Dumper;
use Getopt::Long;
use List::Util qw(sum);

my $varTable;
my $samplesNeeded;
my $sampleControl;
my $cmedianTh = 2;
my $original = 4.3;
my $type = 'somatic';

GetOptions (
           "varTable|v=s"   => \$varTable,
           "samples|s=s"    => \$samplesNeeded,
           "control|c=s"    => \$sampleControl,
           "cmedianTh|m=f"  => \$cmedianTh,
           "original|o=f"   => \$original,
           "type|t=s"       => \$type,
           "help|h"         => sub {
                               print "usage: $0 grep sample table by doing some cleaning work\n\nOptions:\n\t--varTable\tto grep from\n";
                               print "\t--samples\tcomma seperated sample names\n";
                               print "\t--control\tthe name of the control sample, usually just one sample\n";
                               print "\t--cmedianTh\tthe Th for cmedian, default 2\n";
                               print "\t--original\tthe original quality must be greater than, default 4.3\n";
                               print "\t--type\t\teither somatic (default) or germline\n";
                               print "\t--help\t\tprint help\n";
                               print "\n";
                               exit 0;
                             },
           );

my @samples = split(',', $samplesNeeded);
my $samplesRegex = "(";
$samplesRegex .= join "|", @samples;
$samplesRegex .= ")";
print STDERR "sample name regex is $samplesRegex\n";
my $numSamples = scalar(@samples);
print STDERR "number of samples is $numSamples\n";



open IN, "$varTable";
my %colnames;
my %colindex;
my @header;
my $headerPrinted = 0;
my $somindex;
#my $dronindex;
my $gnindex;
my $glindex;
my $gfindex;
my $caddindex;
my $gerpindex;
my $siftindex;
my $polyphenindex;
my $ncindex;
my @ncheader;
my $freqindex;
while ( <IN> ) {
  chomp;
  if (/^[\#]?chr\t/) {
    #it is header
    $_ =~ s/^#//;
    my @cols = split /\t/;
    @header = @cols;
    for(my $i = 0; $i <= $#cols; $i++) {
      $colnames{$cols[$i]} = $i;
      $colindex{$i} = $cols[$i];
    } #header info constructed

    $somindex = ($type eq 'somatic')? $colnames{"somatic"}:$colnames{"trace"};
    #$dronindex = $colnames{"dron"};
    $gnindex = $colnames{"geneName"};
    $glindex = $colnames{"geneLoc"};
    $gfindex = $colnames{"functionalClass"};
    $caddindex = $colnames{"CADD_phred"};
    $gerpindex = $colnames{"GERP_RS"};
    $siftindex = $colnames{"SIFT_score"};
    $polyphenindex = $colnames{"Polyphen2_HVAR_pred"};
    $ncindex = $colnames{$sampleControl."maf"};               #normal maf index
    @ncheader = @header[$ncindex-1..$ncindex+1];
    $freqindex = ($type eq 'somatic')? -1:$colnames{"freq"};   #no freq info for somatic

  } else {

    my @cols = split /\t/;

    if ($cols[$somindex] =~ /$samplesRegex/) {         #originally called

      my @ncinfos = @cols[$ncindex-1..$ncindex+1];
      my @resArraySampleData;
      my @resArraySampleDataHeader;
      foreach my $tsample (@samples) {
        my $tcindex = $colnames{$tsample.'maf'};
        my @arrayTmp = @cols[$tcindex-1..$tcindex+1];
        push(@resArraySampleData, @arrayTmp);
        my @arrayTmp2 = @header[$tcindex-1..$tcindex+1];
        push(@resArraySampleDataHeader,@arrayTmp2);
      }

      my @resVector = ($type eq 'somatic')? (@cols[0..1],@cols[3..5],@resArraySampleData,$cols[$gnindex],$cols[$glindex],$cols[$gfindex],$cols[$caddindex],$cols[$gerpindex],$cols[$siftindex],$cols[$polyphenindex],$cols[$somindex],@ncinfos):(@cols[0..1],@cols[3..5],@resArraySampleData,$cols[$gnindex],$cols[$glindex],$cols[$gfindex],$cols[$caddindex],$cols[$gerpindex],$cols[$siftindex],$cols[$polyphenindex],$cols[$somindex],@ncinfos,$cols[$freqindex]);
      my @resVectorHeader = ($type eq 'somatic')? (@header[0..1],@header[3..5],@resArraySampleDataHeader,$header[$gnindex],$header[$glindex],$header[$gfindex],$header[$caddindex],$header[$gerpindex],$header[$siftindex],$header[$polyphenindex],$header[$somindex],@ncheader):(@header[0..1],@header[3..5],@resArraySampleDataHeader,$header[$gnindex],$header[$glindex],$header[$gfindex],$header[$caddindex],$header[$gerpindex],$header[$siftindex],$header[$polyphenindex],$header[$somindex],@ncheader,$header[$freqindex]);

      my $maxMaf = 0;
      my $maxTqual = 0;
      my $ssb = 0;
      my $ssbc = 0;
      my $ssbp = 0;
      my $refdav = 0;
      my $altdav = 0;
      my $totalMaf = 0;
      my $totalAlt = 0;
      my $totalRef = 0;

      foreach my $tsample ( @samples ) {                 #first round, etermine max maf and qual
        my $tcindex = $colnames{$tsample.'maf'};
        my @ss = split("\\|", $cols[$tcindex]);
        my $mafTmp = $ss[0];
        my $oriTqualTmp = 0;
        if ( $cols[$tcindex-1] =~ /\|/ ) {
          my @ori = split("\\|", $cols[$tcindex-1]);
          if ($type eq 'somatic') {
            my @oriQual = split(",", $ori[2]);
            $oriTqualTmp = $oriQual[0];
          } elsif ($type eq 'germline') {
            $oriTqualTmp = $ori[1];
          }
        }
        if ( $mafTmp > $maxMaf ) {
          $maxMaf = $mafTmp                            #define max maf
        }
        if ( $original > 0 & $oriTqualTmp > $maxTqual ) {
          $maxTqual = $oriTqualTmp
        }
      }                                                #get max maf and qual for all samples

      foreach my $tsample (@samples) {                 #second round, extract info

        my $tcindex = $colnames{$tsample.'maf'};
        my $oriTqual = 0;
        if ( $cols[$tcindex-1] =~ /\|/ ) {
          my @ori = split("\\|", $cols[$tcindex-1]);
          if ( $type eq 'somatic' ) {
            my @oriQual = split(",", $ori[2]);
            $oriTqual = $oriQual[0];
          } elsif ($type eq 'germline') {
            $oriTqual = $ori[1];
          }
        }

        my @ss = split("\\|", $cols[$tcindex]);
        my $mafTmp = $ss[0];
        my $endBias = $ss[1];
        my @cmeme = split(",", $ss[2]);
        my $cmedianav = $cmeme[1];
        my $cmemeSum = sum(@cmeme);

        my @strandBiases = split(",", $ss[3]);
        my $strandBias = $strandBiases[0];
        my $strandBiasRef = 0;
        my $strandBiasFisherP = -1;
        if ($#strandBiases > 0) {
          $strandBiasRef = $strandBiases[1];
          $strandBiasFisherP = $strandBiases[2];
        }
        my $mappingBias = $ss[4];

        #decide a b allele count
        my $refnow = round($cols[$tcindex+1]*(1-$mafTmp));
        my $altnow = round($cols[$tcindex+1]*$mafTmp);

        if ( $mafTmp > 0 ) {
          $ssb = $ssb + $strandBias*$altnow;
          $ssbp = $ssbp + $strandBiasFisherP*$altnow;
          $refdav = $refdav + $refnow;
          $altdav = $altdav + $altnow;
          $ssbc = $ssbc + 1;
        }

        #decide now
        my $mafNow = 0;
        if ( $mafTmp == 0 ) {
          $mafNow = $mafTmp;
        } elsif ($endBias < 0.9 and (($strandBias != 0 and $strandBias != 1) or ($strandBiasFisherP > 0.7 and $refnow >= 10 and $altnow >= 5 and $mafTmp >= 0.1)) and $mappingBias < 0.8 and $cmemeSum < 5.2 and $cmedianav < $cmedianTh) {
          $mafNow = $mafTmp;
          if ($original > 0) {  #if original quality was required
            if ($oriTqual < $original) {
              if (!($maxMaf > 0.1 and $maxTqual >= $original)) { #set for max T qual
                $mafNow = 0;
              }
            }
          }
        } else {                #save some low VAF ones
          if ($maxMaf > 0.2 and ($original == 0 or ($original > 0 and $maxTqual >= $original)) and $mafTmp >= 0.01) {
            if ($cmedianav < 4) {
              $mafNow = $mafTmp;
            } else {
              $mafNow = 0;
            }
          } elsif ($maxMaf > 0.05 and ($original == 0 or ($original > 0 and $maxTqual < $original)) and $mafTmp >= 0.04) {
            if ($cmedianav == 1 and $mappingBias == 0) {
              $mafNow = $mafTmp;
            } else {
              $mafNow = 0;
            }
          } else {
            $mafNow = 0;
          }
        }

        $totalMaf += $mafNow;
        $totalAlt += $altnow;
        $totalRef += $refnow;

        my @resVectorAdd = ($mafNow, 0, 0, 0, 0, $refnow, $altnow);
        push @resVector, @resVectorAdd;
        my @resVectorAddHeader = ($tsample."mafc",$tsample."mafa",$tsample."ccf",$tsample."ccfSD",$tsample."time",$tsample."refc",$tsample."altc");
        push @resVectorHeader, @resVectorAddHeader;

      }                         #second round, extract info

      if ($altdav > 0) {
        $ssb = $ssb/$altdav;
        $ssbp = $ssbp/$altdav;
      }
      if ($ssbc > 0) {
        $altdav = $altdav/$ssbc;
        $refdav = $refdav/$ssbc;
      }

      if ($altdav != 0 and ($ssb > 0.9 or $ssb < 0.1)) {                         #multiple sample strand bias
        next;           #do not print
      }
      if ($totalMaf == 0) {
        next;
      }

      my $mergeMAFC = $totalAlt/($totalAlt+$totalRef);
      my $mergeMAFA = $mergeMAFC;
      my $mergedCCF = $mergeMAFC;
      my $mergedCCFsd = $mergeMAFC;
      push (@resVector, ($totalMaf, $totalAlt, $totalRef, $mergeMAFC, $mergeMAFA, $mergedCCF, $mergedCCFsd));
      push (@resVectorHeader, ("totalMaf", "totalAlt", "totalRef", "mergeMAFC", "mergeMAFA", "mergeCCF", "mergeCCFsd"));

      if ($headerPrinted == 0) {
        printf("%s\n",join("\t", @resVectorHeader));
        $headerPrinted += 1;
      }
      printf("%s\n",join("\t", @resVector));

    } #originally called
  } #other than header line
} #each line

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}


exit 0;
