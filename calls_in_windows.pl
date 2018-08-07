#!/usr/bin/perl
#this is a script for finding the common variations among multiple samples using a sliding window method

use strict;
use Getopt::Long;
use File::Glob ':glob';
use File::Basename;
use Data::Dumper;
use FindBin qw($RealBin);

my $cmdline=join(" ",@ARGV);


my %opt = (
           'chr'          => undef,
           'window'       => undef,
	   'mutation'     => undef,
           'mutationT'    => undef,
           'normal'       => undef,
           'hetero'       => undef,
	   'cd10'         => 1,
           'indel'        => undef,
           'indelT'       => undef,
           'sv'           => undef,
	   'nonrepeat'    => undef,
           'nonrecurrent' => undef,
           'nonselfchain' => undef,
           'nonsnp'       => undef,
	   'denovo'       => undef,
	   'common'       => undef,
           'tmpdir'       => "./",
           'prefix'       => undef,
           'sampleInfo'   => undef,
           'type'         => undef,
           'tolerance'    => 0,
          );

GetOptions (
            "chr|c=s"        => \$opt{'chr'},
            "window|w=i"     => \$opt{'window'},
            "mutation|m=s"   => \$opt{'mutation'},
            "mutationT|x=s"  => \$opt{'mutationT'},
            "normaln=s"      => \$opt{'normal'},
            "hetero|e"       => \$opt{'hetero'},
            "indel|i=s"      => \$opt{'indel'},
            "indelT|y=s"     => \$opt{'indelT'},
            "sv=s"           => \$opt{'sv'},
            "nonrepeat|r=s"  => \$opt{'nonrepeat'},
            "nonrecurrent|u" => \$opt{'nonrecurrent'},
            "nonselfchain|s=s" => \$opt{'nonselfchain'},
            "nonsnp|p"       => \$opt{'nonsnp'},
            "denovo|d"       => \$opt{'denovo'},
            "common|o=s"     => \$opt{'common'},
            "tmpdir=s"       => \$opt{'tmpdir'},
            "prefix=s"       => \$opt{'prefix'},
            "sampleInfo=s"   => \$opt{'sampleInfo'},
            'type=s'         => \$opt{'type'},
            "tolerance=i"    => \$opt{'tolerance'},
            "help|h"         => sub {
              print "usage: $0 [options]\n\nOptions:\n\t--chr\t\tthe chromosome name, like X, 22 etc.. if not set, search for all the chromosomes\n";
              print "\t--window\tthe window size of searching variations\n";
              print "\t--mutation[T]\tonly search for substitutions, 'T' for the already prepared mutation table\n";
              print "\t--normal\tthe sample id for normal or blood samples to decide somatic mutations\n";
              print "\t--hetero\tonly search for hetero substitutions.\n";
              print "\t--indel[T]\tsearch for partial indel and depth indel, 'T' for the already prepared indel table\n";
              print "\t--sv\t\tsearch for SVs\n";
              print "\t--nonrepeat\toptionally choose whether allow variations in repetitive regions to be searched, default, yes, provide rep file\n";
              print "\t--nonselfchain\toptionally choose whether allow variations in selfchain regions to be searched, default, yes, provide sc file\n";
              print "\t--nonrecurrent\tskip calls that are recurrent across samples.\n";
              print "\t--nonsnp\toptionally choose whether allow variations that are known in dbSNP to be searched, default, yes\n";
              print "\t--denovo\twhether substract Parents or not\n";
              print "\t--common\tthe number of the common patients\n";
              print "\t--help\t\tprint this help message\n";
              print "\t--tolerance\tthe distance tolerance to be allowed for comparing indels for generating tables\n";
              print "\t--sampleInfo\tthe sample infomation table (tumor normal each line sep by tab)\n";
              print "\t--type\tthe type of mutation, somatic or germline\n";
              print "\t--tmpdir\tthe temporary dir to write tmp files\n";
              print "\t--prefix\tthe prefix for sample IDs\n";
              print "\n";
              exit 0;
            },
           );


my @prefix = split(',', $opt{'prefix'});
my $prefixReg = join('|', @prefix);
print STDERR "prefixReg is $prefixReg\n";

print STDERR "gathering mutation type is $opt{type}\n";

my %somatic;
my %germline;   #may have multiple tumors
if ($opt{'sampleInfo'}) {
  open IN, "$opt{sampleInfo}" or die "$opt{sampleInfo} is not readable.\n";
  while ( <IN> ) {
    chomp;
    s/[\s\n]$//;
    my @columns = split /\t/;
    my $tumor = $columns[0];
    my $normal = $columns[1];
    $somatic{$tumor} = $normal;
    push(@{$germline{$normal}}, $tumor) if $normal ne 'undef';
  }
  close IN;
  print STDERR "somatic Info is:\n";
  print STDERR Dumper (\%somatic);
  #print STDERR Dumper (\%germline);
}


my $bin = $RealBin;
print "@ $0 $cmdline\n" if $opt{'window'};

my %normals;
foreach my $normalSample (split(',', $opt{'normal'})){
  $normals{$normalSample} = '';
}
print STDERR "normal samples:\n";
print STDERR Dumper(\%normals);

my %variations;
my $winsize = $opt{window};
my @common;
while ($opt{common} =~ /(\d+)/g) {
  push @common, $1;
}


my %iupac;
$iupac{R}{A}="G";
$iupac{R}{G}="A";
$iupac{Y}{C}="T";
$iupac{Y}{T}="C";
$iupac{S}{C}="G";
$iupac{S}{G}="C";
$iupac{W}{A}="T";
$iupac{W}{T}="A";
$iupac{M}{A}="C";
$iupac{M}{C}="A";
$iupac{K}{G}="T";
$iupac{K}{T}="G";


#############################repeat region loaded in 2D hash %repeatmask############################################
my %repeatmask;
my @rs_rm;                      #start array for each chr
my $old_chr_rm;                 #checking the chr
my $ptr_rm;                     #pointer for repeatmask sub
my $old_run_rm;

if ($opt{nonrepeat}) {
  open REP, "$opt{nonrepeat}";
  while (<REP>) {
    next if /^#/;
    chomp;
    my @tmp = split /\t/;
    (my $chr = $tmp[0]) =~ s/^chr(\w+)$/$1/;
    my $repeat_s = $tmp[3];
    my $repeat_e = $tmp[4];
    $repeatmask{$chr}{$repeat_s} = $repeat_e;
  }
  close REP;
  print STDERR "#repeatmasker file loaded\n";
}


#############################Ucsc Self Chain %selfChain#############################################
my %selfChain;
my @rs_selfChain;
my $old_chr_selfChain;
my $ptr_selfChain;
my $old_run_selfChain;

if ($opt{nonselfchain}) {
  open SELFCHAIN, "$opt{nonselfchain}";
  while ( <SELFCHAIN> ) {
    next if /^#/;
    chomp;
    my ($bin, $score, $tName, $tSize, $tStart, $tEnd, $qName, $qSize, $qStrand, $qStart, $qEnd, $id, $normScore) = split /\t/;

    $tName =~ s/^chr(\w+)$/$1/;
    $qName =~ s/^chr(\w+)$/$1/;

    next if (($normScore > 0 and $normScore < 20) or $score < 2000);

    next if ($selfChain{$tName}{$tStart} ne '' and $selfChain{$tName}{$tStart} >= $tEnd);
    $selfChain{$tName}{$tStart} = $tEnd;

    next if ($selfChain{$qName}{$qStart} ne '' and $selfChain{$qName}{$qStart} >= $qEnd);
    $selfChain{$qName}{$qStart} = $qEnd;
  }
  close SELFCHAIN;
  print STDERR "#selfchain annotation loaded\n";
}


#############################common snp region loaded in 2D hash %snpmask############################################
my %snpmask;
my @rs_snp;      #start array for each chr
my $old_chr_snp; #checking the chr
my $ptr_snp;     #pointer for repeatmask sub
my $old_run_snp;

if ($opt{nonsnp}) {
  open SNP, "";
  while (<SNP>) {
    next if /^#/;
    chomp;
    my @tmp = split /\t/;
    (my $chr = $tmp[0]) =~ s/^chr(\w+)$/$1/;
    my $snp_s = $tmp[3];
    my $snp_e = $tmp[4];
    $snpmask{$chr}{$snp_s} = $snp_e;
  }
  close SNP;
  print STDERR "#snp130 file loaded\n";
}


####################################load the file of substitution with the high density score#########################################################
my %snv;
my %samples;
my %netcoh;
if ($opt{mutation}) {
  open MUT, "$opt{mutation}";
  my @snv_files;
  while (<MUT>) {
    chomp;
    push(@snv_files, $_);
  }
  close MUT;

  foreach my $snv_file (@snv_files) {
    my $individual;

    $snv_file =~ /($opt{'prefix'}\d+)[^0-9a-zA-Z]/;
    $individual = $1;

    open SNV, "$snv_file";
    my $revertornot = "no";
    my $printerror = 0;
    my $singlecalling = "no";
    while ( <SNV> ) {
      chomp;
      if ($_ =~ /^#/) {
        if ($_ =~ /^#CHROM\tPOS\tID/) {
          my @cols = split (/\t/, $_);
          my $minusI = 1;
          print STDERR "$cols[$#cols]\n";
          if ( exists($normals{$cols[$#cols - $minusI]}) ) {
            $revertornot = "yes";
          } elsif ( $cols[$#cols - $minusI] eq 'FORMAT' ) {
            $singlecalling = "yes";
          }
          print STDERR "revert or not: $revertornot\n";
          print STDERR "singlecalling: $singlecalling\n";
          next;
        } else {
          next;
        }
      }

      my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $blood, $sample) = split /\t/;

      if ($revertornot eq 'yes') {   #revert sample and blood
        my $tmp = $sample;
        $sample = $blood;
        $blood = $tmp;
      }

      my @formats = split(':', $FORMAT);
      my %formindex;
      for(my $f = 0; $f <= $#formats; $f++) {
        $formindex{$formats[$f]} = $f;
      }
      if ($printerror == 0) {
        print STDERR Dumper(\%formindex);
        $printerror ++;
      }

      $CHROM =~ s/chr//;
      my $function;
      my $MAF;
      if ( !$opt{chr} or ($opt{chr} and ($CHROM eq $opt{chr})) ) {

        if ($QUAL < 30){
          next;                          #skip low quality calls
        }

        if ($ID ne '.' or $INFO =~ /dbSNP/ or $INFO =~ /1KG/ or $INFO =~ /ESP5400\=/) {  #snp in population, is it a somatic one?

          my $freq = -1;
          if ($INFO =~ /(1KG=(.+?));/) {
            my $kid = $1;       #re define $id when absent
            $freq = $2;
            if ($ID eq '.') {
              $ID = $kid;
            }
          }

          if ($INFO =~ /(ESP5400=(.+?));/) {
            my $kid = $1;       #re define $id to ESP when absent
            $freq = $2;
            if ($ID eq '.') {
              $ID = $kid;
            }
          }

          #judging somatic or not
          my $somatic = 0;
          if ($blood ne '') {
            my @blood = split(/\:/, $blood);
            my @bad;
            if (exists($formindex{'AD'})) {   #AD found
              @bad = split (/\,/, $blood[$formindex{'AD'}]);
            }
            if ($blood[$formindex{'GT'}] !~ /1/ and $bad[1] == 0) {
              $somatic = 1;
            }
          }

          if ($somatic == 0) { #keep somatic ones even if it is marked as a common snp
            if ($freq == -1) {
              next;
            } elsif ($freq > 0 and ($ID !~ /^1KG/ and $ID !~ /^ESP5400/)) {
              next;
            }
          }
        }                       #somatic common snp

      PRODUCE:

        if ($INFO =~ /MQ0Fraction=(.+?);/) {
          next if ($1 > 0.1);   #skip multiple matching region
        }
        if ($INFO =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/) {
          my $strandb = $1;
          my $baseqb = $2;
          my $mapqb = $3;
          my $tailb = $4;
          if ($strandb =~ /e/) { #strand bias
            next;
          } elsif ($strandb < 0.005) { #strand bias
            next;
          } elsif ($baseqb =~ /e/) { #basequality bias
            next;
          } elsif ($baseqb < 0.0005) { #basequality bias
            next;
          } elsif ($mapqb =~ /e/) { #mapquality bias
            next;
          } elsif ($mapqb < 0.0001) { #mapquality bias
            next;
          } elsif ($tailb =~ /e/) { #tailbias
            next;
          } elsif ($tailb < 0.005) { #tailbias
            next;
          } else {
            #pass
          }
        }

        my $coor = $CHROM.':'.$POS;

        $INFO =~ /[\;]?DP\=(\d+)[\;]?/;
        my $dp = $1;
        next if ($dp < 10);              #skip low coverage region
        $INFO =~ /(function=.+?)$/;
        $function = $1;
        $INFO =~ /[\;]?DP4\=((\d+)\,(\d+)\,(\d+)\,(\d+))[\;]?/;
        my $depth_record = $dp.'('.$1.')';
        my $depth_var = $4+$5;
        next if $depth_var < 2;
        $MAF = sprintf("%.3f", $depth_var/($2+$3+$4+$5));
        $sample =~ /^([01]\/[01])\:/;
        if ($1 eq '0/1' or $1 eq '1/0') {
          $variations{$CHROM}{$POS}{$individual}{'SUB'}{'info'} = 'SUB_Hetero:'.$CHROM.':'.$POS.':'.$REF.'->'.$ALT.':'.$depth_record.':'.$function;
        } elsif ($1 eq '1/1') {
          $variations{$CHROM}{$POS}{$individual}{'SUB'}{'info'} = 'SUB_Homo:'.$CHROM.':'.$POS.':'.$REF.'->'.$ALT.':'.$depth_record.':'.$function;
        }

        if ($opt{'table'}) {   #generate the mutation table
          $snv{$coor}{$individual} = $MAF;
          $snv{$coor}{'function'} = $function;
          $snv{$coor}{'info'} = join("\t", ($ID,$REF,$ALT));
        }
      } #the defined chromosome
    } #each SNV file
    close SNV;
    print STDERR "#$snv_file loaded\n";
    $samples{$individual} = "";
    $netcoh{$individual} = "";
  }
  print STDERR "#all SNVs' files loaded\n";
}


if ($opt{mutationT}) {
  my %header;
  open MUT, "$opt{mutationT}";
  while ( <MUT> ){
    chomp;
    my @cols;
    if ($_ =~ /^[^#]?chr\t/){  #header line
      @cols = split /\t/;
      for (my $i = 0; $i <= $#cols; $i++){
        $header{$i} = $cols[$i];
      }
    } #remember header
    else {
      @cols = split /\t/;
      my $chr;
      my $pos;
      my $id;
      my $ref;
      my $alt;
      my $function;
      my $rep;
      my $sc;
      my %sampleNeed;
      if ($opt{type}){
        for (my $i = 0; $i <= $#cols; $i++) {
          if ($header{$i} eq $opt{type}) {
            next if ($cols[$i] eq 'NA');
            my @somaticSamps = split(',', $cols[$i]);
            foreach my $somaticSamp (@somaticSamps) {
              $somaticSamp =~ s/\[\w+\]//;
              next unless (exists $somatic{$somaticSamp});
              $sampleNeed{$somaticSamp}++;
            }
          }
        }
      } #judge which sample is needed

      for (my $i = 0; $i <= $#cols; $i++){
        if ($header{$i} eq 'chr'){
          $chr = $cols[$i];
        }
        elsif ($header{$i} eq 'pos'){
          $pos = $cols[$i];
        }
        elsif ($header{$i} eq 'id'){
          $id = $cols[$i];
        }
        elsif ($header{$i} eq 'ref'){
          $ref = $cols[$i];
        }
        elsif ($header{$i} eq 'alt'){
          $alt = $cols[$i];
        }
        elsif ($header{$i} =~ /function/) {
          $function = $cols[$i];
        }
        elsif ($header{$i} eq 'rep') {
          $rep = $cols[$i];
        }
        elsif ($header{$i} eq 'sc') {
          $sc = $cols[$i];
        }
        elsif ($header{$i} =~ /^($prefixReg)(\d+)?([A-Za-z0-9\-\_]+)?maf/) {
          my $individual = $header{$i};
          (my $individualWithoutMAF = $individual) =~ s/maf$//;
          next if (exists $normals{$individualWithoutMAF});
          next if (!exists $somatic{$individualWithoutMAF});
          if ($opt{type}){
            next if (!exists $sampleNeed{$individualWithoutMAF});
          }
          my $maf = $cols[$i];
          if ($maf =~ /^([0-9\.]+)\|/){
            next if $1 <= 0.04;  #skip where maf less than 0.04
          } else {
            next if $maf == 0;
          }
          my $depth = $cols[$i+1];
          $samples{$individualWithoutMAF} = "";
          $netcoh{$individualWithoutMAF} = "";
          next if ($opt{'nonsnp'} and $id != ".");
          $variations{$chr}{$pos}{$individualWithoutMAF}{'SUB'}{'info'} = 'SNV:'.$chr.':'.$pos.':'.$id.':'.$ref.'->'.$alt.':'.$maf.':'.$depth.':'.$function;
          $variations{$chr}{$pos}{'rep'} = $rep;
          $variations{$chr}{$pos}{'sc'} = $sc;
        }
      }
    }
  }
  close MUT;
}


####################################load the file of indel###############################################################
my %indel;
if ($opt{indel}) {
  open INDEL, "$opt{indel}";
  my @indel_files;
  while ( <INDEL> ) {
    chomp;
    push(@indel_files, $_);
  }
  close INDEL;

  foreach my $indel_file (@indel_files) {
    my $individual;
    if ($indel_file =~ /ARJ/) {
      $individual = "AC3";
      if ($indel_file =~ /type1/){
         $individual .= "T1";
      } elsif ($indel_file =~ /type2/) {
         $individual .= "T2";
      } else {
         $individual .= "U";
      }
      open ARJ, "$indel_file";
      while ( <ARJ> ) {
        chomp;
        next if /^chr\t/;
        next if /^#/;
        my @cols = split /\t/;
        my $CHROM = $cols[0];
        my $POS = $cols[1];
        my $coor = $CHROM.':'.$POS;
        my $MAF = 0;

        if ($indel_file =~ /type1/){
          for(my $i = 9; $i <= 42; $i = $i+3) {
            $MAF += $cols[$i];
          }
          $MAF = sprintf("%.3f", $MAF/12);
        } elsif ($indel_file =~ /type2/){
          $MAF = sprintf("%.3f", ($cols[45]+$cols[48])/2);
        } else {
          my $startI = ($opt{'clinical'})? 5:9;
          my $endI = ($opt{'clinical'})? 18:48;
          my $increI = ($opt{'clinical'})? 1:3;
          for(my $i = $startI; $i <= $endI; $i = $i+$increI) {
            $MAF += $cols[$i];
          }
          $MAF = sprintf("%.3f", $MAF/14);
        }

        my $function = $cols[5];
        if ($opt{'clinical'}){
          $function = $cols[19];
        }
        $variations{$CHROM}{$POS}{$individual}{'INDEL'}{'info'} = $CHROM.':'.$POS.'(maf='.$MAF.';'.$function.')';
        $variations{$CHROM}{$POS}{$individual}{'INDEL'}{'end'} = $POS + 1;
        if ($opt{'table'}) {   #generate the mutation table
           $indel{$coor}{$individual} = $MAF;
           $indel{$coor}{'function'} = $function;
           $indel{$coor}{'info'} = join("\t", ($cols[2], $cols[3], $cols[4]));
           if ($opt{'clinical'}){
             $snv{$coor}{'clinical'} = $cols[20];
           }
        }
      }
      close ARJ;
      print STDERR "#ARJ $individual loaded.\n";
      $samples{$individual} = "";
      next;
    } #ARJ files

    $indel_file =~ /($opt{'prefix'}\d+)[^0-9a-zA-Z]/;
    $individual = $1;

    my %recheck;
    my $recheckN = 0;
    if ($opt{'recheck'}) {   #do recheck here
      my @rechecks = bsd_glob("$opt{recheck}/indel/compared_*/$individual");
      $recheckN = scalar(@rechecks);
      print STDERR "$recheckN\n";
      foreach my $recheck (@rechecks) {
        open RECHECK, "$recheck";
        while ( <RECHECK> ) {
           chomp;
           my @cols = split /\t/;
           my $coor = $cols[0].':'.$cols[1];
           if ($cols[6] == 0 and $cols[5] > 7) {                       #only record for absent indels
              $recheck{$coor} += 1;
            }
         }
        close RECHECK;
      } #each recheck file
    } #do recheck

    if ($opt{'clinical'}) { # want to grep the clinical sites only, first generate tmp intersect file
      my $cmd = "perl $bin/intersectFiles.pl -o $indel_file -m $opt{'clinical'} -vcf -overlap -column INFO -t $opt{'tolerance'} >$opt{'tmpdir'}/tmp";
      RunCommand($cmd,0,0);
      $indel_file = "$opt{'tmpdir'}/tmp";
    }

    open COHORT, "$indel_file";
    my $revertornot = "no";
    my $printerror = 0;
    my $singlecalling = "no";
    while ( <COHORT> ) {
      chomp;
      if ($_ =~ /^#/) {
        if ($_ =~ /^#CHROM\tPOS\tID/) {
          my @cols = split /\t/;
          my $minusI = ($opt{'clinical'})? 2:1;
          print STDERR "$cols[$#cols - 1]\n";
          if ( exists($normals{$cols[$#cols - $minusI]}) ) {
            $revertornot = "yes";
          } elsif ($cols[$#cols - $minusI] eq 'FORMAT') {
            $singlecalling = "yes";
          }
          print STDERR "revert or not: $revertornot\n";
          print STDERR "singlecalling: $singlecalling\n";
          next;
        } else {
          next;
        }
      }

      my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $sample, $blood, $clinINFO) = split /\t/;

      if ($revertornot eq 'yes') {   #revert sample and blood
        my $tmp = $sample;
        $sample = $blood;
        $blood = $tmp;
      }

      my @formats = split(':', $FORMAT);
      my %formindex;
      for(my $f = 0; $f <= $#formats; $f++) {
        $formindex{$formats[$f]} = $f;
      }
      if ($printerror == 0){
        print STDERR Dumper(\%formindex);
        $printerror ++;
      }


      $CHROM =~ s/chr//;
      my $function;
      my $MAF;
      if ( !$opt{chr} or ($opt{chr} and ($CHROM eq $opt{chr})) ) {

        if ($QUAL < 30) {
          next;                          #skip low quality calls
        }

         if ($opt{'clinical'}){  # want to grep the clinical sites only
           goto PRODUCEINDEL;
         }

        if ($ID ne '.' or $INFO =~ /dbSNP/ or $INFO =~ /1KG/ or $INFO =~ /ESP5400\=/) {  #indel in population, is it a somatic one?

          my $freq = -1;
          if ($INFO =~ /(1KG=(.+?));/) {
            my $kid = $1;       #re define $id when absent
            $freq = $2;
            if ($ID eq '.') {
              $ID = $kid;
            }
          }

          if ($INFO =~ /(ESP5400=(.+?));/) {
            my $kid = $1;       #re define $id to ESP when absent
            $freq = $2;
            if ($ID eq '.') {
              $ID = $kid;
            }
          }

          #judging somatic or not
          my $somatic = 0;
          if ($blood ne '') {
            my @blood = split(/\:/,$blood);
            my @bad;
            if (exists($formindex{'AD'})) {   #AD found
              @bad = split (/\,/, $blood[$formindex{'AD'}]);
            }

            if (scalar(@bad) > 0){
              if ($blood[$formindex{'GT'}] !~ /1/ and $bad[1] == 0) {
                $somatic = 1;
              }
            } else {
              if ($blood[$formindex{'GT'}] !~ /1/) {
                $somatic = 1;
              }
            }
          }

          if ($somatic == 0) { #keep somatic ones even if it is marked as a common snp
            if ($freq == -1) {
              next;
            } elsif ($freq > 0 and ($ID !~ /^1KG/ and $ID !~ /^ESP5400/)) {
              next;
            }
          }
        }                       #somatic common snp

      PRODUCEINDEL:

        if ($INFO =~ /MQ0Fraction=(.+?);/) {
          next if ($1 > 0.1);            #skip multiple matching region
        }
        if ($INFO =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/) {
          my $strandb = $1;
          my $baseqb = $2;
          my $mapqb = $3;
          my $tailb = $4;
          if ($strandb =~ /e/) { #strand bias
            next;
          } elsif ($strandb < 0.005) { #strand bias
            next;
          } elsif ($baseqb =~ /e/) { #basequality bias
            next;
          } elsif ($baseqb < 0.0005) { #basequality bias
            next;
          } elsif ($mapqb =~ /e/) { #mapquality bias
            next;
          } elsif ($mapqb < 0.0001) { #mapquality bias
            next;
          } elsif ($tailb =~ /e/) { #tailbias
            next;
          } elsif ($tailb < 0.005) { #tailbias
            next;
          } else {
            #pass
          }
        }


        my $coor = $CHROM.':'.$POS;
        next if ($opt{'recheck'} and (!exists($recheck{$coor}) or $recheck{$coor} < $recheckN));   #skip the rechecked ones


        $INFO =~ /[\;]?DP\=(\d+)[\;]?/;
        my $dp = $1;
        next if ($dp < 10);          #skip low coverage region
        $INFO =~ /(function=.+?)$/;
        $function = $1;
        $INFO =~ /[\;]?DP4\=((\d+)\,(\d+)\,(\d+)\,(\d+))[\;]?/;
        my $depth_record = $dp.'('.$1.')';
        my $depth_var = $4+$5;
        next if $depth_var < 2;
        $MAF = sprintf("%.3f", $depth_var/($2+$3+$4+$5));
        $variations{$CHROM}{$POS}{$individual}{'INDEL'}{'info'} = 'INDEL:'.$CHROM.':'.$POS.':'.$REF.'->'.$ALT.':'.$depth_record.':'.$function;
        $variations{$CHROM}{$POS}{$individual}{'INDEL'}{'end'} = $POS + 1;

        if ($opt{'table'}) {   #generate the mutation table
          $indel{$coor}{$individual} = $MAF;
          $indel{$coor}{'function'} = $function;
          $indel{$coor}{'info'} = join("\t", ($ID, $REF, $ALT));
          if ($opt{'clinical'}) {
            $indel{$coor}{'clinical'} = $clinINFO; #clinical id information
          }
        }  #generate table

      } #the same chromosome
    } #each cohort
    close COHORT;
    print STDERR "#$indel_file loaded\n";
    $samples{$individual} = "";
    $netcoh{$individual} = "";
  }
  print STDERR "#all indels' files loaded\n";
}



if ($opt{indelT}) {
  my %header;
  open IND, "$opt{indelT}";
  while ( <IND> ){
    chomp;
    my @cols;
    if ($_ =~ /^[#]?chr\t/){  #header line
      @cols = split /\t/;
      for (my $i = 0; $i <= $#cols; $i++){
        $header{$i} = $cols[$i];
      }
    } #remember header
    else {
      @cols = split /\t/;
      my $chr;
      my $pos;
      my $id;
      my $ref;
      my $alt;
      my $function;
      my $rep;
      my $sc;
      my %sampleNeed;
      if ($opt{type}){
        for (my $i = 0; $i <= $#cols; $i++) {
          if ($header{$i} eq $opt{type}) {
            next if ($cols[$i] eq 'NA');
            my @somaticSamps = split(',', $cols[$i]);
            foreach my $somaticSamp (@somaticSamps) {
              $somaticSamp =~ s/\[\w+\]//;
              next unless (exists $somatic{$somaticSamp});
              $sampleNeed{$somaticSamp}++;
            }
          }
        }
      } #judge which sample is needed

      for (my $i = 0; $i <= $#cols; $i++){
        if ($header{$i} eq 'chr'){
          $chr = $cols[$i];
        }
        elsif ($header{$i} eq 'pos'){
          $pos = $cols[$i];
        }
        elsif ($header{$i} eq 'id'){
          $id = $cols[$i];
        }
        elsif ($header{$i} eq 'ref'){
          $ref = $cols[$i];
        }
        elsif ($header{$i} eq 'alt'){
          $alt = $cols[$i];
        }
        elsif ($header{$i} eq 'func'){
          $function = $cols[$i];
        }
        elsif ($header{$i} eq 'rep') {
          $rep = $cols[$i];
        }
        elsif ($header{$i} eq 'sc') {
          $sc = $cols[$i];
        }
        elsif ($header{$i} =~ /^($prefixReg)(\d+)?([A-Za-z0-9\-\_]+)?maf/){
          my $individual = $header{$i};
          (my $individualWithoutMAF = $individual) =~ s/maf$//;
          next if ($normals{$individualWithoutMAF});
          next if (!exists $somatic{$individualWithoutMAF});
          if ($opt{type}){
            next if (!exists $sampleNeed{$individualWithoutMAF});
          }
          my $maf = $cols[$i];
          if ($maf =~ /^([0-9\.]+)\|/){
            next if $1 <= 0.04;  #skip where maf less than 0.04
          } else {
            next if $maf == 0;
          }
          my $depth = $cols[$i+1];
          $samples{$individualWithoutMAF} = "";
          $netcoh{$individualWithoutMAF} = "";
          next if ($opt{'nonsnp'} and $id != ".");
          $variations{$chr}{$pos}{$individualWithoutMAF}{'INDEL'}{'info'} = 'INDEL:'.$chr.':'.$pos.':'.$id.':'.$ref.'->'.$alt.':'.$maf.':'.$depth.':'.$function;
          $variations{$chr}{$pos}{'rep'} = $rep;
          $variations{$chr}{$pos}{'sc'} = $sc;
        }
      }
    }
  }
  close IND;
}


####################################load the file of SV###############################################################
my %sv;
if ($opt{'sv'}) {
  open SV, "$opt{'sv'}";
  my @sv_files;
  while ( <SV> ) {
    chomp;
    push(@sv_files, $_);
  }
  close SV;

  foreach my $sv_file (@sv_files) {

    my $individual;

    $sv_file =~ /($opt{'prefix'}\d+)[^0-9a-zA-Z]/;
    $individual = $1;

    open COHORT, "$sv_file";
    my $printerror = 0;
    while ( <COHORT> ) {
      chomp;
      if ($_ =~ /^#/) {
          next;
      }

      my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $blood, $sample) = split /\t/;

      my @formats = split(':', $FORMAT);
      my %formindex;
      for(my $f = 0; $f <= $#formats; $f++) {
        $formindex{$formats[$f]} = $f;
      }
      if ($printerror == 0){
        print STDERR Dumper(\%formindex);
        $printerror ++;
      }

      $CHROM =~ s/chr//;
      my $function;
      my $MAF;
      if ( !$opt{chr} or ($opt{chr} and ($CHROM eq $opt{chr})) ) {

        my $coor = $CHROM.':'.$POS;

        $INFO =~ /[\;]?SR\=(\d+)[\;]?/;
        my $sr = $1;
        $INFO =~ /[\;]?RP\=(\d+)[\;]?/;
        my $rp = $1;
        $INFO =~ /[\;]?REF\=(\d+)[\;]?/;
        my $refsu = $1;
        $INFO =~ /[\;]?REFPAIR\=(\d+)[\;]?/;
        my $refpair = $1;

        my $vard = $sr + $rp;
        my $depth = $sr + $rp + $refsu + $refpair;
        $MAF = sprintf("%.3f", $vard/$depth);
        $variations{$CHROM}{$POS}{$individual}{'SV'}{'info'} = 'SV:'.$CHROM.':'.$POS.':'.$REF.'->'.$ALT.':'.$MAF.':('.$sr.'+'.$rp.'+'.$refsu.'+'.$refpair.')';
        $variations{$CHROM}{$POS}{$individual}{'SV'}{'end'} = $POS + 1;

      } #the same chromosome
    } #each cohort
    close COHORT;
    print STDERR "#$sv_file loaded\n";
    $samples{$individual} = "";
  }
  print STDERR "#all SV files loaded\n";
}

print STDERR Dumper(\%samples);


#################################generate a big array with all the starts############################################################
foreach my $chr (sort keys %variations) { #foreach chromosome $chr

  my @pre_starts = sort {$a<=>$b} keys %{$variations{$chr}};
  my @starts;

  foreach my $start (@pre_starts) {

    my @n_individuals = keys %{$variations{$chr}{$start}};
    my $n_individuals = scalar(@n_individuals);

    if ($opt{nonrecurrent}) {
      if ($n_individuals >= 2) {                                       #skip recurrent vars
         next;
      }
    }
    if ($opt{nonselfchain}) {
      if ($opt{'mutationT'} or $opt{'indelT'}){
        next if ($variations{$chr}{$start}{'sc'} == 1);
      } else {
        next if (selfChainMask('run1', $chr, $start) == 1);
      }
    }

    if ($opt{nonrepeat}) {
      if ($opt{'mutationT'} or $opt{'indelT'}) {
        next if ($variations{$chr}{$start}{'rep'} == 1);
      } else {
        next if (repeatmask('run1', $chr, $start, $start) == 1);
      }
    }

    push @starts, $start;

  }

  my $totalwinnum = scalar @starts;
  print STDERR "#total window chr$chr: $totalwinnum\n";


  ###################################window starts to move!!!!! go go go ######################################################

  my $ptr = 0;       #pointer for moving the index of the starts array
  my $old_start_ptr = -1; #pointer of the start coor of the last window
  my %old_starts;         #checking for the old starts;
  my %result;
  my $n_result_windows = 0;
  my $n_result_regions = 0;
  my $last_found_window_end = 0;
  my $last_region_start = 0;
  my %merged = ();

  foreach my $start (@starts) {

    my $winend = $start+$winsize-1;
    my $winregion = 'chr'.$chr.':'.$start."-".$winend; #decide the window covered region

    #delete the old variation when the window moving one by one
    unless ( $old_start_ptr < 0 ) {

      #get all the old starts retained in the result
      foreach my $os (sort {$a<=>$b} keys %old_starts) {
        my $flag_for_delete_os = 1;
        foreach my $individual (keys %{$variations{$chr}{$os}}) {
          foreach my $var_type (keys %{$variations{$chr}{$os}{$individual}}) {
            if ($var_type =~ /SUB/) {
              delete $result{$individual}{$os}{$var_type};
            } else { #the variation is an indel, then check whether the indel overlaps the current window
              my $end = $variations{$chr}{$os}{$individual}{$var_type}{end};
              if ($end < $start) { #nonoverlapping, delete it
                delete $result{$individual}{$os}{$var_type};
              }
            }
          }
          if (keys %{$result{$individual}{$os}} == 0) {
            delete $result{$individual}{$os};
          } else {          #if there is any one still have this start
            $flag_for_delete_os = 0;
          }

          if (keys %{$result{$individual}} == 0) {
            delete $result{$individual};
          }
        }
        delete $old_starts{$os} if ($flag_for_delete_os == 1);

      }  #foreach
    }    #unless

    #which variations should be added into this window?
    while ($ptr <= $#starts and $starts[$ptr] <= $winend) { #now scan the starts array within the window region

      #generate the result hash whose keys is the individual names
      my $start_here = $starts[$ptr];
      foreach my $individual (keys %{$variations{$chr}{$start_here}}) {
        foreach my $var_type (keys %{$variations{$chr}{$start_here}{$individual}}) {
          $result{$individual}{$start_here}{$var_type} = $variations{$chr}{$start_here}{$individual}{$var_type}{'info'};
        }
      }
      $ptr++;
    }

    $old_start_ptr++; #move one by one of the old start ptr for deleting the keys
    $old_starts{$starts[$old_start_ptr]} = ''; #add the old starts;

    my $individuals = join "..", (sort keys %result);

    my $np = 0;
    while ($individuals =~ /$prefixReg/g) {
      $np++;
    }
    print STDERR "$chr\t$start\t$np\n";

    #checking and output

    if ((scalar(@common)==1 and $np >= $common[0]) or (scalar(@common)==2 and $np >= $common[0] and $np <= $common[1])) {

      if ($start > $last_found_window_end) { #this window is NOT overlapping with previous window
        $n_result_regions++;  #add region
        $merged{$start}{'end'} = $winend;
        foreach my $individual_here (sort keys %result) {
          foreach my $start_here (sort keys %{$result{$individual_here}}) {
            foreach my $var_type_here (sort keys %{$result{$individual_here}{$start_here}}) {
              $merged{$start}{$individual_here}{$result{$individual_here}{$start_here}{$var_type_here}} = '';
            }
          }
        } #add to merge
        $last_region_start = $start;
      } else { # an old region, add window
        $merged{$last_region_start}{'end'} = $winend;
        foreach my $individual_here (sort keys %result) {
          foreach my $start_here (sort keys %{$result{$individual_here}}) {
            foreach my $var_type_here (sort keys %{$result{$individual_here}{$start_here}}) {
              $merged{$last_region_start}{$individual_here}{$result{$individual_here}{$start_here}{$var_type_here}} = '';
            }
          }
        } #add to merge
      } # an old region
      $n_result_windows++;
      $last_found_window_end = $winend;
    }

    my $window_number = $old_start_ptr+1;
  }
  print STDERR "#Number of total_variation_starts/nonoverlapping_regions come up chr$chr: $n_result_windows/$n_result_regions\n";

  #print
  foreach my $region_s (sort {$a <=> $b} keys %merged) {
     my $region_e = $merged{$region_s}{'end'};
     my $region_size = $region_e - $region_s + 1;
     print "chr$chr\:$region_s\-$region_e\t$region_size\n";
     foreach my $individual (sort keys %{$merged{$region_s}}) { #individual
        next if $individual eq 'end';
        my $info = "$individual\t";
        foreach my $variant (sort keys %{$merged{$region_s}{$individual}}){
          $info .= $variant.';';
        }
        print "$info\n";
     } #individual
  }

}

exit 0;

##########################################################################################################
#                                        Subroutine Region                                               #
##########################################################################################################
sub repeatmask {
    my ($run, $chr, $start, $end) = @_;
    my $flag = 0;
    if (($chr ne $old_chr_rm) or ($run ne $old_run_rm)){
       @rs_rm = sort {$a <=> $b} keys %{$repeatmask{$chr}};
       $ptr_rm = 0;
    }
    while (($ptr_rm<=$#rs_rm) and ($repeatmask{$chr}{$rs_rm[$ptr_rm]} < $start)){
      $ptr_rm++;
    }
    if ($rs_rm[$ptr_rm] <= $end){
      $flag = 1;
    }
    $old_chr_rm = $chr;
    $old_run_rm = $run;
    return $flag;
}

sub snpmask {
    my ($run, $chr, $start, $end) = @_;
    my $flag = 0;
    if (($chr ne $old_chr_snp) or ($run ne $old_run_snp)){
       @rs_snp = sort {$a <=> $b} keys %{$snpmask{$chr}};
       $ptr_snp = 0;
    }
    while (($ptr_snp<=$#rs_snp) and ($snpmask{$chr}{$rs_snp[$ptr_snp]} < $start)){
      $ptr_snp++;
    }
    if ($rs_snp[$ptr_snp] <= $end){
      $flag = 1;
    }
    $old_chr_snp = $chr;
    $old_run_snp = $run;
    return $flag;
}

sub selfChainMask {
    my ($run, $chr, $coor) = @_;
    my $flag = 0;
    if (($chr ne $old_chr_selfChain) or ($run ne $old_run_selfChain)){
        @rs_selfChain = sort {$a <=> $b} keys %{$selfChain{$chr}};
        $ptr_selfChain = 0;
    }
    while (($ptr_selfChain <= $#rs_selfChain) and ($selfChain{$chr}{$rs_selfChain[$ptr_selfChain]} < $coor)){
        $ptr_selfChain++;
    }
    if ($rs_selfChain[$ptr_selfChain] <= $coor){
        $flag = 1;
    }
    $old_chr_selfChain = $chr;
    $old_run_selfChain = $run;
    return $flag;
}


sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}
