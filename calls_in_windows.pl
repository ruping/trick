#!/usr/bin/perl
#this is a script for finding the common variations among multiple samples using a sliding window method

use strict;
use Getopt::Long;
use File::Glob ':glob';
use File::Basename;
use Data::Dumper;

my $cmdline=join(" ",@ARGV);
print "@ $0 $cmdline\n";


my %opt = (
           'chr'          => undef,
           'window'       => undef,
	   'mutation'     => undef,
           'mutationT'    => undef,
           'hetero'       => undef,
	   'cd10'         => 1,
           'indel'        => undef,
           'indelT'       => undef,
	   'nonrepeat'    => undef,
           'nonrecurrent' => undef,
           'nonselfchain' => undef,
           'nonsnp'       => undef,
	   'denovo'       => undef,
	   'common'       => undef,
           'table'        => undef,
           'recheck'      => undef,
          );

GetOptions (
           "chr|c=s"        => \$opt{chr},
           "window|w=i"     => \$opt{window},
	   "mutation|m=s"   => \$opt{mutation},
           "mutationT|x=s"  => \$opt{mutationT},
           "hetero|e"       => \$opt{hetero},
           "indel|i=s"      => \$opt{indel},
           "indelT|y=s"     => \$opt{indelT},
	   "nonrepeat|r"    => \$opt{nonrepeat},
           "nonrecurrent|u" => \$opt{nonrecurrent},
           "nonselfchain|s" => \$opt{nonselfchain},
           "nonsnp|p"       => \$opt{nonsnp},
	   "denovo|d"       => \$opt{denovo},
	   "common|o=s"     => \$opt{common},
           "table|t"        => \$opt{table},
           "recheck|k=s"    => \$opt{recheck},
	   "help|h"         => sub{
	                       print "usage: $0 [options]\n\nOptions:\n\t--chr\t\tthe chromosome name, like X, 22 etc.. if not set, search for all the chromosomes\n";
                               print "\t--window\tthe window size of searching variations\n";
                               print "\t--mutation[T]\tonly search for substitutions, 'T' for the already prepared mutation table\n";
			       print "\t--hetero\tonly search for hetero substitutions.\n";
                               print "\t--indel[T]\t\tsearch for partial indel and depth indel, 'T' for the already prepared indel table\n";
                               print "\t--nonrepeat\toptionally choose whether allow variations in repetitive regions to be searched, default, yes\n";
                               print "\t--nonselfchain\toptionally choose whether allow variations in selfchain regions to be searched, default, yes\n";
                               print "\t--nonrecurrent\tskip calls that are recurrent across samples.\n";
                               print "\t--nonsnp\toptionally choose whether allow variations that are known in dbSNP to be searched, default, yes\n";
			       print "\t--denovo\twhether substract Parents or not\n";
			       print "\t--common\tthe number of the common patients\n";
			       print "\t--help\t\tprint this help message\n";
                               print "\t--table\t\tgenerate table of variants\n";
                               print "\t--recheck\tthe dir whether recheck files are located under ./compared_\*/\n";
                               print "\n";
                               exit 0;
       	                     },
           );


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
  open REP, "/ifs/data/c2b2/ac_lab/rs3412/annotation/hg19/hg19.repeats_UCSC.gff";
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
  open SELFCHAIN, "/ifs/data/c2b2/ac_lab/rs3412/annotation/hg19/hg19.SelfChain_UCSC.txt";
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
    if ($snv_file =~ /ARJ/) {
      $individual = "AC3";
      if ($snv_file =~ /type1/){
         $individual .= "T1";
      } elsif ($snv_file =~ /type2/) {
         $individual .= "T2";
      } else {
         $individual .= "U";
      }
      open ARJ, "$snv_file";
      while ( <ARJ> ) {
        chomp;
        next if /^AC/;
        my @cols = split /\t/;
        my $CHROM = $cols[0];
        my $POS = $cols[1];
        my $coor = $CHROM.':'.$POS;
        my $MAF = 0;
        if ($snv_file =~ /type1/){
          for(my $i = 9; $i <= 42; $i = $i+3) {
            $MAF += $cols[$i-1];
          }
          $MAF = sprintf("%.3f", $MAF/12);
        } elsif ($snv_file =~ /type2/) {
          $MAF = sprintf("%.3f", ($cols[44]+$cols[47])/2);
        } else {
          for(my $i = 9; $i <= 48; $i = $i+3){
            $MAF += $cols[$i-1];
          }
          $MAF = sprintf("%.3f", $MAF/14);
        }
        my $function = $cols[49];
        $variations{$CHROM}{$POS}{$individual}{'SUB'}{'info'} = $CHROM.':'.$POS.'(maf='.$MAF.';'.$function.')';
        if ($opt{'table'}) {   #generate the mutation table
           $snv{$coor}{$individual} = $MAF;
           $snv{$coor}{'function'} = $function;
           $snv{$coor}{'info'} = join("\t", ($cols[2], $cols[3], $cols[4]));
        }
      }
      close ARJ;
      print STDERR "#ARJ $individual loaded.\n";
      $samples{$individual} = "";
      next;
    }
    $snv_file =~ /(AC\d+)[^0-9a-zA-Z]/;
    $individual = $1;

    my %recheck;
    my $recheckN = 0;
    if ($opt{'recheck'}){   #do recheck here
      my @rechecks = bsd_glob("$opt{recheck}/snv/compared_*/$individual");
      $recheckN = scalar(@rechecks);
      print STDERR "$recheckN\n";
      foreach my $recheck (@rechecks) {
        open RECHECK, "$recheck";
        while ( <RECHECK> ) {
           chomp;
           my @cols = split /\t/;
           my $coor = $cols[0].':'.$cols[1];
           if (($cols[4]+$cols[5]+$cols[6]+$cols[7]) == 0 and $cols[2] > 7) {                       #only record for absent indels
              $recheck{$coor} += 1;
           }
        }
        close RECHECK;
      } #each recheck file
    } #do recheck

    open SNV, "$snv_file";
    while ( <SNV> ) {
      next if /^#/;
      chomp;
      my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $GT1, $GT2) = split /\t/;

      $CHROM =~ s/chr//;
      my $function;
      my $MAF;
      if ( !$opt{chr} or ($opt{chr} and ($CHROM eq $opt{chr})) ) {

        if ($QUAL < 30){
          next;                          #skip low quality calls
        }

        if ($ID ne '.' or $INFO =~ /dbSNP/ or $INFO =~ /1KG/) {  #snp in population, is it a somatic one?
          my $somatic = 0;

          my $freq = -1;
          if ($INFO =~ /(1KG=(.+?));/) {
            my $kid = $1;       #re define $id when absent
            $freq = $2;
            if ($ID eq '.') {
              $ID = $kid;
            }
          }

          if ($GT2 ne ''){
            my @formats = split(/\:/, $FORMAT);
            my $ADindex = -1;
            for (my $i = 0; $i <= $#formats; $i++){
              if ($formats[$i] eq 'AD'){
                $ADindex = $i;
                last;
              }
            }
            my @blood = split(/\:/,$GT2);
            my @bad;
            if ($ADindex != -1) {   #AD found
              @bad = split (/\,/, $blood[$ADindex]);
            }
            if ($blood[0] eq '0/0' and $bad[1] == 0) {
              $somatic = 1;
            }
          }

          if ($somatic == 0) { #keep somatic ones even if it is marked as a common snp
            if ($freq == -1) {
              next;
            } elsif ($freq > 0.0005) {
              next;
            }
          }
        }                       #somatic common snp


        if ($INFO =~ /MQ0Fraction=(.+?);/) {
          next if ($1 > 0.1);            #skip multiple matching region
        }
        if ($INFO =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/) {
          my $tailb = $4;
          if ($tailb =~ /e/) {
            next;                        #skip tail bias region
          } elsif ($tailb < 0.005){
            next;
          }
        }


        my $coor = $CHROM.':'.$POS;
        next if ($opt{'recheck'} and (!exists($recheck{$coor}) or $recheck{$coor} < $recheckN));

        $INFO =~ /\;DP\=(\d+);/;
        my $dp = $1;
        next if ($dp < 10);              #skip low coverage region
        $INFO =~ /(function=.+?)$/;
        $function = $1;
        $INFO =~ /\;DP4\=((\d+)\,(\d+)\,(\d+)\,(\d+))\;/;
        my $depth_record = $dp.'('.$1.')';
        my $depth_var = $4+$5;
        next if $depth_var < 2;
        $MAF = sprintf("%.3f", $depth_var/($2+$3+$4+$5));
        $GT1 =~ /^([01]\/[01])\:/;
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
    if ($_ =~ /^chr\t/){  #header line
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
        elsif ($header{$i} =~ /^AC\d+maf$/ || $header{$i} =~ /^AC3[TU]/){
          my $individual = $header{$i};
          next if ($individual eq 'AC1maf' or $individual eq 'AC3maf' or $individual eq 'AC547maf' or $individual eq 'AC581maf');
          #next unless ($individual =~ /AC3[TU]/ or $individual eq 'AC546maf' or $individual eq 'AC580maf');
          my $maf = $cols[$i];
          next if $maf < 0.05;  #skip where maf less than 0.05
          my $depth;
          if ($individual eq 'AC3T1'){$depth = $cols[$i+4]; $samples{$individual} = "";}
          elsif ($individual eq 'AC3T2'){$depth = $cols[$i+3]; $samples{$individual} = "";}
          elsif ($individual eq 'AC3U'){$depth = $cols[$i+2]; $samples{$individual} = "";}
          else{$depth = $cols[$i+1]; $samples{$individual} = ""; $netcoh{$individual} = "";}
          $variations{$chr}{$pos}{$individual}{'SUB'}{'info'} = 'SNV:'.$chr.':'.$pos.':'.$id.':'.$ref.'->'.$alt.':'.$maf.':'.$depth.':'.$function;
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
        next if /^chr/;
        my @cols = split /\t/;
        my $CHROM = $cols[0];
        my $POS = $cols[1];
        my $coor = $CHROM.':'.$POS;
        my $MAF = 0;

        if ($indel_file =~ /type1/){
          for(my $i = 9; $i <= 42; $i = $i+3) {
            $MAF += $cols[$i-1];
          }
          $MAF = sprintf("%.3f", $MAF/12);
        } elsif ($indel_file =~ /type2/){
          $MAF = sprintf("%.3f", ($cols[44]+$cols[47])/2);
        } else {
          for(my $i = 9; $i <= 48; $i = $i+3) {
            $MAF += $cols[$i-1];
          }
          $MAF = sprintf("%.3f", $MAF/14);
        }

        my $function = $cols[49];
        $variations{$CHROM}{$POS}{$individual}{'INDEL'}{'info'} = $CHROM.':'.$POS.'(maf='.$MAF.';'.$function.')';
        $variations{$CHROM}{$POS}{$individual}{'INDEL'}{'end'} = $POS + 1;
        if ($opt{'table'}) {   #generate the mutation table
           $indel{$coor}{$individual} = $MAF;
           $indel{$coor}{'function'} = $function;
           $indel{$coor}{'info'} = join("\t", ($cols[2], $cols[3], $cols[4]));
        }
      }
      close ARJ;
      print STDERR "#ARJ $individual loaded.\n";
      $samples{$individual} = "";
      next;
    }
    $indel_file =~ /(AC\d+)[^0-9a-zA-Z]/;
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

    open COHORT, "$indel_file";
    while ( <COHORT> ) {
      next if /^#/;
      chomp;
      my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $GT1, $GT2) = split /\t/;

      $CHROM =~ s/chr//;
      my $function;
      my $MAF;
      if ( !$opt{chr} or ($opt{chr} and ($CHROM eq $opt{chr})) ) {

        if ($QUAL < 30) {
          next;                          #skip low quality calls
        }

        if ($ID ne '.' or $INFO =~ /dbSNP/ or $INFO =~ /1KG/) {  #snp in population, is it a somatic one?
          my $somatic = 0;

          my $freq = -1;
          if ($INFO =~ /(1KG=(.+?));/) {
            my $kid = $1;       #re define $id when absent
            $freq = $2;
            if ($ID eq '.') {
              $ID = $kid;
            }
          }

          if ($GT2 ne '') {
            my @formats = split(/\:/, $FORMAT);
            my $ADindex = -1;
            for (my $i = 0; $i <= $#formats; $i++){
              if ($formats[$i] eq 'AD'){
                $ADindex = $i;
                last;
              }
            }
            my @blood = split(/\:/,$GT2);
            my @bad;
            if ($ADindex != -1) {   #AD found
              @bad = split (/\,/, $blood[$ADindex]);
            }
            if (scalar(@bad) > 0){
              if ($blood[0] eq '0/0' and $bad[1] == 0) {
                $somatic = 1;
              }
            } else {
              if ($blood[0] eq '0/0') {
                $somatic = 1;
              }
            }
          }

          if ($somatic == 0) { #keep somatic ones even if it is marked as a common snp
            if ($freq == -1) {
              next;
            } elsif ($freq > 0.0005) {
              next;
            }
          }
        }                       #somatic common snp

        if ($INFO =~ /MQ0Fraction=(.+?);/) {
          next if ($1 > 0.1);            #skip multiple matching region
        }
        if ($INFO =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/) {
          my $tailb = $4;
          if ($tailb =~ /e/) {
            next;                        #skip tail bias region
          } elsif ($tailb < 0.005){
            next;
          }
        }


        my $coor = $CHROM.':'.$POS;
        next if ($opt{'recheck'} and (!exists($recheck{$coor}) or $recheck{$coor} < $recheckN));   #skip the rechecked ones


        $INFO =~ /\;DP\=(\d+);/;
        my $dp = $1;
        next if ($dp < 10);          #skip low coverage region
        $INFO =~ /(function=.+?)$/;
        $function = $1;
        $INFO =~ /\;DP4\=((\d+)\,(\d+)\,(\d+)\,(\d+))\;/;
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
    if ($_ =~ /^chr\t/){  #header line
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
        elsif ($header{$i} =~ /^AC\d+maf$/ || $header{$i} =~ /^AC3[TU]/){
          my $individual = $header{$i};
          next if ($individual eq 'AC1maf' or $individual eq 'AC3maf' or $individual eq 'AC547maf' or $individual eq 'AC581maf');
          #next unless ($individual =~ /AC3[TU]/ or $individual eq 'AC546maf' or $individual eq 'AC580maf');
          my $maf = $cols[$i];
          next if $maf <= 0.05;  #skip where maf less than 0.05
          my $depth;
          if ($individual eq 'AC3T1'){$depth = $cols[$i+4]; $samples{$individual} = "";}
          elsif ($individual eq 'AC3T2'){$depth = $cols[$i+3]; $samples{$individual} = "";}
          elsif ($individual eq 'AC3U'){$depth = $cols[$i+2]; $samples{$individual} = "";}
          else{$depth = $cols[$i+1]; $samples{$individual} = ""; $netcoh{$individual} = "";}
          $variations{$chr}{$pos}{$individual}{'INDEL'}{'info'} = 'INDEL:'.$chr.':'.$pos.':'.$id.':'.$ref.'->'.$alt.':'.$maf.':'.$depth.':'.$function;
        }
      }
    }
  }
  close IND;
}





##generate big table

if ($opt{'table'} and $opt{'mutation'}){
  print "#chr\tpos\tid\tref\talt";
  foreach my $name (sort keys %samples) {
    print "\t$name";
  }
  print "\tfunction\n";
}
if ($opt{'table'} and $opt{'indel'}){
  print "#chr\tpos\tid\tref\talt";
  foreach my $name (sort keys %samples) {
    print "\t$name";
  }
  print "\tfunction\n";
}



#################################generate a big array with all the starts############################################################
foreach my $chr (sort keys %variations) { #foreach chromosome $chr

  my @pre_starts = sort {$a<=>$b} keys %{$variations{$chr}};
  my @starts;

  foreach my $start (@pre_starts) {
    my @n_individuals = keys %{$variations{$chr}{$start}};
    my $n_individuals = scalar(@n_individuals);

    if ($opt{nonrecurrent}) {
      if ($n_individuals >= 2) {                                       #skip recurrent point mutation
         my $recuflag = 0;
         foreach my $indi (@n_individuals){
            $recuflag = 1 if exists($netcoh{$indi});
         }
         next if $recuflag == 1;
      }
    }
    if ($opt{nonselfchain}) {
      next if (selfChainMask('run1', $chr, $start) == 1);
    }

    if ($opt{nonrepeat} and $opt{nonsnp}) {
      unless (repeatmask('run1', $chr, $start, $start) == 1 or snpmask('run1', $chr, $start, $start) == 1) {
        push @starts, $start;
      }
    } elsif ($opt{nonrepeat} and !$opt{nonsnp}) {
      unless (repeatmask('run1', $chr, $start, $start) == 1){
        push @starts, $start;
      }
    } elsif (!$opt{nonrepeat} and !$opt{nonsnp}) {
      push @starts, $start;
    }
  }


  #print table for snv
  if ($opt{'table'} and $opt{'mutation'}) {
    foreach my $start (@starts) {
      my $coor = $chr.":".$start;
      my $info = $snv{$coor}{'info'};
      print "$chr\t$start\t$info";
      foreach my $name (sort keys %samples) {
        if ($snv{$coor}{$name} ne '') {
          print "\t$snv{$coor}{$name}";
        } else {
          print "\t0";
        }
      }
      my $function = $snv{$coor}{'function'};
      print "\t$function\n";
    }
    goto TABLE;
  }
  #print table for indel
  if  ($opt{'table'} and $opt{'indel'}){
    foreach my $start (@starts) {
      my $coor = $chr.":".$start;
      my $info = $indel{$coor}{'info'};
      print "$chr\t$start\t$info";
      foreach my $name (sort keys %samples) {
        if ($indel{$coor}{$name} ne '') {
          print "\t$indel{$coor}{$name}";
        } else {
          print "\t0";
        }
      }
      my $function = $indel{$coor}{'function'};
      print "\t$function\n";
    }
    goto TABLE;
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
          $result{$individual}{$start_here}{$var_type} = $variations{$chr}{$start_here}{$individual}{$var_type}{info};
        }
      }
      $ptr++;
    }

    $old_start_ptr++; #move one by one of the old start ptr for deleting the keys
    $old_starts{$starts[$old_start_ptr]} = ''; #add the old starts;

    my $individuals = join "..", (sort keys %result);

    my $np = 0;
    while ($individuals =~ /AC/g) {
      $np++;
    }
    print STDERR "$chr\t$start\t$np\n";

    #checking and output

    if ((scalar(@common)==1 and $np >= $common[0]) or (scalar(@common)==2 and $np >= $common[0] and $np <= $common[1])) {
      #my $flag = 0;
      #foreach my $individual_here ( sort keys %result ) {
      #  my $info;
      #  foreach my $start_here (sort keys %{ $result{$individual_here}} ) {
      #    foreach my $var_type_here (sort keys %{$result{$individual_here}{$start_here}}) {
      #      $info .= $result{$individual_here}{$start_here}{$var_type_here}.';';
      #    }
      #  }
      #  print "\t\t\t$individual_here\t$info\n" if ($flag == 1);
      #  if ($flag == 0) {
      #    print "$winregion\t$individual_here\t$info\n";
      #    $flag++;
      #  }
      #}                         #foreach individual
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

TABLE:

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
