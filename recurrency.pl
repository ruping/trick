use strict;
use Data::Dumper;
use File::Glob ':glob';
use List::Util qw[min max];
use Data::Dumper;
use Getopt::Long;

my $file;
my $task;
my $type;
my $somaticInfo;        # if for somatic judgement
my $bloodCall;          # whether the blood is also single called
my $rna;                # if for rna found

my $Th_tumorLOD = 5.3;
my $Th_normalLOD = 2.3;
my $errorRate = 0.001;
my $Th_maf = 0.03;
my $Th_endsratio = 0.9;
my $Th_vard = 2;
my $Th_badQualFrac = 0.6;
my $Th_cmeancmedian = 5.5;
my $Th_cmedian = 2;
if ($type eq 'indel') {
  $Th_endsratio = 0.95;
  $Th_cmeancmedian = 6.5;
  $Th_cmedian = 3;
}


GetOptions (
           "file|f=s"          => \$file,             #filename
           "type|t=s"          => \$type,             #snv or indel
           "task|k=s"          => \$task,             #task type
           "somaticInfo|s=s"   => \$somaticInfo,      #info for somatic sample pairs
           "bloodCall|b=s"     => \$bloodCall,        #bloodCall
           "Th_cmeme=f"        => \$Th_cmeancmedian,  #cmeancmedian
           "rna|r=s"           => \$rna,              #rna info
           "help|h"         => sub {
                               print "usage: $0 final preparation of somatic/germline variant calls \n\nOptions:\n\t--file\t\tthe filename of the res file\n";
                               print "\t--type\t\tthe type of variants, snv or indel\n";
                               print "\t--task\t\tthe aim of the analysis, maf, somatic, filter, (sam)founds, (rna)trace, depth, freq, etc\n";
                               print "\t--somaticInfo\tthe paired sample information\n";
                               print "\t--bloodCall\twhether the blood has called for variants (yes or no), for a sanity check\n";
                               print "\t--rna\t\tthe dna-rna matching table, when task is rnatrace\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );


print STDERR "Th_maf: $Th_maf\n";
print STDERR "Th_endsratio: $Th_endsratio\n";
print STDERR "Th_vard: $Th_vard\n";
print STDERR "Th_badQualFrac: $Th_badQualFrac\n";
print STDERR "Th_cemancmedian: $Th_cmeancmedian\n";
print STDERR "Th_cmedian: $Th_cmedian\n";


my @all;
my %somatic;
my %germline;  #may have multiple tumors
if ($somaticInfo and -s "$somaticInfo") {

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
  #print STDERR Dumper (\%somatic);
  #print STDERR Dumper (\%germline);

}


my %rnasamps;
if ($rna) {
  if ( (! -e $rna) and $rna =~ /,/ ) {
    my @rnas = split(',', $rna);
    foreach my $rnas (@rnas) {
      $rnasamps{$rnas} = '';
    }
  } elsif (-e $rna) {
    my %rcol;
    open IN, "$rna";
    while ( <IN> ) {
      chomp;
      my @cols = split /\t/;
      if ($_ =~ /^([\#])?dna\t/) { #header
        for (my $i=0; $i<=$#cols; $i++) {
          $rcol{$cols[$i]} = $i;
        }
      } else {
        next if $cols[$rcol{'rna'}] eq 'NA';
        $rnasamps{$cols[$rcol{'rna'}]} = $cols[$rcol{'dna'}];
      }
    }
    close IN;
  }
  print STDERR Dumper(\%rnasamps);
}


open IN, "$file";
my %colnames;
my %colindex;
while ( <IN> ) {
  chomp;
  if (/^[\#]?chr\t/) {
    #it is header
    $_ =~ s/^#//;
    my @cols = split /\t/;
    for(my $i = 0; $i <= $#cols; $i++) {
      $colnames{$cols[$i]} = $i;
      $colindex{$i} = $cols[$i];
    }
    if ($task eq 'maf') {
      print "$_\tmaf\n";
    } elsif ($task eq 'freq') {
      print "$_\tfreq\n";
    } elsif ($task eq 'trace') {
      print "$_\ttrace\n";
    } elsif ($task eq 'rnatrace') {
      print "$_\trnatrace\n";
    } elsif ($task eq 'founds') {
      print "$_\tfounds\n";
    } elsif ($task eq 'rnafounds') {
      print "$_\trnafounds\n";
    } elsif ($task =~ /filter/) {
      print "$_\tfilter\n";
    } elsif ($task =~ /somatic/) {
      print "$_\tsomatic\tgermline\n";
    } elsif ($task eq 'depth'){
      print "$_\tdepthav\n";
    } elsif ($task eq 'samfounds'){
      print "$_\tfounds\trsam\n";
    } else {
      print STDERR "task if wierd\n";
      exit 22;
    }
  } else {
    my @cols = split /\t/;
    if ($task eq 'maf') {    #get real maf
      my $maf = 0;
      my $sampleCounts = scalar(@all);
      foreach my $sample (@all) {
        if ($cols[$colnames{$sample}] >= 0.1) {
          my $vard = sprintf("%.1f", $cols[$colnames{$sample}]*$cols[$colnames{$sample}+1]);
          if ($vard >= 2) {
            $maf += $cols[$colnames{$sample}];
          }
        }
      }
      $maf = sprintf("%.6f",$maf/$sampleCounts);
      print "$_\t$maf\n";
    } elsif ($task eq 'freq') {
      my $freq = 'NA';
      my $function = (exists($colnames{'function'}))? $cols[$colnames{'function'}] : die("no function column.\n");
      my $id = $cols[$colnames{'id'}];
      if ($id =~ /1KG=(.+?)$/){
        $freq = $1;
      } elsif ($id =~ /ESP\d+=(.+?)$/){
        $freq = $1;
      } elsif ($function =~ /1KG=(.+?)\;/){
        $freq = $1;
      } elsif ($function =~ /ESP\d+=(.+?)\;/){
        $freq = $1;
      }
      print "$_\t$freq\n";
    } elsif ($task =~ /founds/ or $task =~ /trace/) {    #trace all samples and check whether it is originally called
      my $founds = 0;
      my $samfounds = 0;
      my $trace = '';
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colindex{$i} =~ /^(.+?)maf$/) {
          my $samp = $1;
          if ($task =~ /rna/) {      #for rna
            next if (!exists($rnasamps{$samp}));
          }
          my $maf = $cols[$i];
          my $endsratio = 0;
          my $cmean = 0;
          my $cmedian = 0;
          my $strandRatio = 0;
          my $badQualFrac = 0;

          if ($cols[$i] =~ /\|/) { #split the var surrounding information
            my @infos = split(/\|/, $cols[$i]);
            $maf = $infos[0];
            $endsratio = $infos[1];
            ($cmean, $cmedian) = split(',', $infos[2]);
            $strandRatio = $infos[3];
            $badQualFrac = $infos[4];
          }

          my $depth = $cols[$i+1];
          if ($depth =~ /\,/){
            my @depths = split(',', $cols[$i+1]);   #spandepth, jumpdepth
            $depth = $depths[0];
          }
          my $vard = sprintf("%.1f", $maf*$depth);

          if (($endsratio <= $Th_endsratio or ((1-$endsratio)*$vard >= $Th_vard)) and $badQualFrac <= $Th_badQualFrac and ($strandRatio > 0 and $strandRatio < 1) and (($cmean+$cmedian) < $Th_cmeancmedian or $cmedian <= $Th_cmedian)) {  #it looks good
            if ($task =~ /rna/) {
              if ($maf >= ($Th_maf - 0.02) and $vard >= $Th_vard) {
                $founds++;
                $trace .= "$samp,";
              }
            } else {
              if ($maf >= $Th_maf and $vard >= ($Th_vard+1)) {
                if ($somaticInfo ne '') { #count only tumor
                  if ( exists($somatic{$samp}) ) {
                    $founds++;
                    $trace .= "$samp,";
                    if ($colindex{$i-1} =~ /^$samp$/) { #calling information for the same sample
                      if ( $cols[$i-1] =~ /\|/ ) { #it is called originally
                        $samfounds++;
                      }
                    }
                  }
                } else {
                  $founds++;                       #count all samples
                  $trace .= "$samp,";
                  if ($colindex{$i-1} =~ /^$samp$/) { #calling information for the same sample
                    if ( $cols[$i-1] =~ /\|/ ) { #it is called originally
                      $samfounds++;
                    }
                  }
                }
              }                 #maf and vard requirements
            } #for dna
          }  # it looks good
        } #maf
      } #each column
      if ($task =~ /founds/) {
        print "$_\t$founds";
        if ($task =~ /sam/) {    #original founds
          my $rsam = ($founds > 0)? sprintf("%.2f", $samfounds/$founds) : 0;
          print "\t$rsam";
        }
        print "\n";
      } elsif ($task =~ /trace/) {
        if ($trace ne ''){
          print "$_\t$trace\n";
        } else {
          print "$_\tNA\n";
        }
      }
    } elsif ($task eq 'depth') {   #find av depth
      my $dep = 0;
      my $Ndep = 0;
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colindex{$i} =~ /maf$/) {
          if ($cols[$i] >= 0.1){    #found clonal
            my $vard = sprintf("%.1f", $cols[$i]*$cols[$i+1]);
            if ($vard >= $Th_vard) {
              $Ndep++;
              $dep += $cols[$i+1];
            }
          }
        } #maf
      } #each column
      my $depav = ($Ndep > 0)? sprintf("%.1f", $dep/$Ndep):0;
      print "$_\t$depav\n";
    } elsif ($task =~ /filter/) {     #filter (based on Jie's Results)
      my @detectedSample;
      my %detectedSample;
      my $somaticCalled = 0;
      my %mafs;
      my $chr;
      my $pos;
      my $endsratio = 0;
      my $strandRatio = 0;
      my $badQualFrac = 0;
      my $cmean = 0;
      my $cmedian = 0;
      my $cmeanav = 0;
      my $cmedianav = 0;
      my $mmaf = 0;
      my $rep = 0;
      my $sc = 0;
      for ( my $i = 0; $i <= $#cols; $i++ ) {
        if ($colindex{$i} eq 'chr') {
          $chr = $cols[$i];
        } elsif ($colindex{$i} eq 'pos') {
          $pos = $cols[$i];
        } elsif ($colindex{$i} eq 'trace') {
          my ($traceSomatic, $traceGermline) = split(';', $cols[$i]);
          $traceSomatic =~ /somatic\=(.+?)$/;
          $traceSomatic = $1;
          $traceGermline =~ /germline\=(.+?)$/;
          $traceGermline = $1;
          if ($traceSomatic =~ /\,/) {
            my @traceSomatic = split(/\,/, $traceSomatic);
            push(@detectedSample, @traceSomatic);
            $somaticCalled = 1;
          }
          if ($traceGermline =~ /\,/) {
            my @traceGermline = split(/\,/, $traceGermline);
            push(@detectedSample, @traceGermline);
          }
          foreach my $detectedSamp (@detectedSample) {
            $detectedSample{$detectedSamp} = '';
          }
          #print STDERR "$chr\t$pos\t$traceSomatic\t$traceGermline\n";
          #print STDERR Dumper(\%detectedSample);
        } elsif ($colindex{$i} eq 'rep') {
          $rep = $cols[$i];
        } elsif ($colindex{$i} eq 'sc') {
          $sc = $cols[$i];
        } elsif ($colindex{$i} =~ /^(.+?)maf$/) {   #store maf information
          my $samp = $1;
          $mafs{$samp} = $cols[$i];
        } elsif ($colindex{$i} eq 'cmeanav') {
          $cmeanav = $cols[$i];
          $cmedianav = $cols[$i+1];
        }
      } #each column

      foreach my $samp (keys %mafs) {  #get endsratio and cmean cmedian info
        next if (! exists($detectedSample{$samp}));   #only look at the sample where it is called
        my $sampmaf = $mafs{$samp};
        if ($sampmaf =~ /\|/) { #split the var surrounding information
          my @infos = split(/\|/, $sampmaf);
          $endsratio = ($infos[0] > $mmaf)? $infos[1]:$endsratio;
          ($cmean, $cmedian) = ($infos[0] > $mmaf)? split(',', $infos[2]):($cmean, $cmedian);
          $strandRatio = ($infos[0] > $mmaf)? $infos[3]:$strandRatio;
          $badQualFrac = ($infos[0] > $mmaf)? $infos[4]:$badQualFrac;
          $mmaf = ($infos[0] > $mmaf)? $infos[0]:$mmaf;
          #print STDERR "sr: $strandRatio\n";
        }
      }

      my $status;
      #print STDERR "$chr\t$pos\t$rep$sc\t$detectedSample[0]\t$mmaf\t$endsratio\t$cmean\t$cmedian\t$cmeanav\t$cmedianav\n";
      if ($rep == 1 and $sc == 1) {
        $status = ($endsratio < $Th_endsratio and $badQualFrac <= ($Th_badQualFrac-0.1) and (($cmean+$cmedian) < ($Th_cmeancmedian-1) or $cmedian < $Th_cmedian) and ($cmeanav + $cmedianav) <= ($Th_cmeancmedian-0.3))? 'PASS':'FOUT';   #conservative for rep and sc
      } elsif ($rep == 1 or $sc == 1) {
        $status = ($endsratio < $Th_endsratio and $badQualFrac <= $Th_badQualFrac and (($cmean+$cmedian) <= $Th_cmeancmedian or $cmedian <= $Th_cmedian) and ($cmeanav + $cmedianav) <= $Th_cmeancmedian)? 'PASS':'FOUT';
      } else {
        $status = ($endsratio < $Th_endsratio and $badQualFrac <= $Th_badQualFrac and (($cmean+$cmedian) <= $Th_cmeancmedian or $cmedian <= $Th_cmedian) and ($cmeanav + $cmedianav) <= $Th_cmeancmedian)? 'PASS':'FOUT';
      }

      if ($type eq 'indel' and $somaticCalled == 1){
        $status = 'PASS';
      }
      if ($task !~ /keep/) {
        print "$_\t$status\n" if ($status eq 'PASS');
        print STDERR "$_\t$status\n" if ($status eq 'FOUT');
      } else {
        print "$_\t$status\n";
      }
    } elsif ($task =~ /somatic/) {  #find somatic ones
      my %tumor;
      my %blood;
      my %bloodcalled;   #store the called information for blood
      my %nonblood;
      my %unknown;
      for ( my $i = 0; $i <= $#cols; $i++ ) {

        if ($colindex{$i} =~ /^(.+?)maf$/) {
          my $samp = $1;
          my $maf = $cols[$i];
          my $endsratio = 0;
          my $strandRatio = 0;
          my $badQualFrac = 0;
          my $cmean = 0;
          my $cmedian = 0;

          if ($cols[$i] =~ /\|/) {  #split the var surrounding information
            my @infos = split(/\|/, $cols[$i]);
            $maf = $infos[0];
            $endsratio = $infos[1];
            ($cmean, $cmedian) = split(',', $infos[2]);
            $strandRatio = $infos[3];
            $badQualFrac = $infos[4];
          }

          my $depth = $cols[$i+1];
          my $vard = sprintf("%.1f", $maf*$depth);

          if (exists $somatic{$samp}) {     #for tumor samples require some additional thing
            if (($endsratio <= $Th_endsratio or ((1-$endsratio)*$vard >= $Th_vard)) and $badQualFrac <= $Th_badQualFrac and ($strandRatio > 0 and $strandRatio < 1) and (($cmean+$cmedian) < ($Th_cmeancmedian-0.3) or $cmedian <= $Th_cmedian)) { #true event
              if ( $maf >= 0.1 ) {  #clonal ones
                $maf = $maf;
              } else {              #subclonal ones, subject to additional constrains
                if ($badQualFrac < ($Th_badQualFrac-0.2) and $cmedian <= ($Th_cmedian-0.3)) {
                  $maf = $maf;
                } else {
                  $maf = 0;         #not reliable somatic
                }
              }
            } else {
              $maf = 0;         #not reliable somatic
            }
          }

          #print STDERR "$samp\t$maf\t$endsratio\t$cmean\t$cmedian\n";

          if (exists $somatic{$samp}) {       #it is tumor sample name
            if ($vard >= ($Th_vard+1) and $maf >= ($Th_maf+0.01) and $depth >= 8) {
              if ($vard >= 50) {              #must be sig anyway
                $tumor{$samp} = $maf;
              } else {                        #cal for less vard
                print STDERR "$samp\t$errorRate\t$maf\t$vard\t$depth\t";
                ##############################################################################
                my $tumorLOD = &calTumorLOD($errorRate, $maf, $vard, $depth); #cal tumor LOD
                ##############################################################################
                print STDERR "$tumorLOD\n";
                if ($tumorLOD >= $Th_tumorLOD) {
                  $tumor{$samp} = $maf;
                }
              }
            }
          } elsif (exists $germline{$samp}) { #it is blood
            foreach my $ct (@{$germline{$samp}}) {
              if ($bloodCall eq 'yes' and $cols[$i-1] =~ /\|/) { #it is originally called
                $bloodcalled{$ct} = '';
              }
              print STDERR "$samp\t$errorRate\t$maf\t$vard\t$depth\t";
              ##############################################################################
              my $normalLOD = &calNormalLOD($errorRate, $maf, $vard, $depth); #cal tumor LOD
              ##############################################################################
              print STDERR "$normalLOD\n";
              if ( $maf == 0 ) {   #no blood alt found
                if ( $normalLOD > $Th_normalLOD and $depth >= 8 ) {  #around 8 depth
                  $nonblood{$ct} = '';
                } else {
                  $unknown{$ct} = '';
                }
              } else {          #maf != 0
                $blood{$ct} = $maf.','.$normalLOD;
              }
            }
          }                     #it is blood
        }                       #maf
      }                         #each column

      my $soma = 'NA';
      my $germ = 'NA';
      foreach my $tumorSamp (keys %tumor) {   #check blood to confirm that it is somatic
        my $stype = 'NA';
        if (exists $nonblood{$tumorSamp}) {
          $stype = 'good';
          $stype .= ($tumor{$tumorSamp} < 0.1)? 'Sub' : '';  #add subclonal info
          $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';

        } elsif (exists($blood{$tumorSamp})) {

          my @bloodmafLOD = split(',',$blood{$tumorSamp});
          if ( $bloodmafLOD[0] <= 0.05 and $tumor{$tumorSamp}/$bloodmafLOD[0] >= 4 and $bloodmafLOD[1] > $Th_normalLOD ) {
            $stype = 'doubt';
            $stype .= ($tumor{$tumorSamp} < 0.1)? 'Sub' : '';  #add subclonal info
            $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';
          } elsif (exists($unknown{$tumorSamp})) {
            $stype = 'undef';
            $stype .= ($tumor{$tumorSamp} < 0.1)? 'Sub' : '';  #add subclonal info
            $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';
          } else {
            $germ = ($germ eq 'NA')? $tumorSamp.',':$germ.$tumorSamp.',';
          }

        } else { #unknown ones, all germ
          $germ = ($germ eq 'NA')? $tumorSamp.',':$germ.$tumorSamp.',';
        }
      }

      if ($soma eq 'NA' and $germ eq 'NA') {     #check whether it is germline
        foreach my $ct (sort keys %bloodcalled){
          $germ .= $ct.',' if ($ct ne '');
        }
      }
      print "$_\t$soma\t$germ\n" if ($soma ne 'NA' or $germ ne 'NA');
    } #somatic
  }
}
close IN;


#calculate LOD scores for deciding somatic calls
sub calTumorLOD {
  my $e = shift;
  my $f = shift;
  my $v = shift;
  my $d = shift;
  my $r = $d-$v;
  my $Pr = $f*($e/3) + (1-$f)*(1-$e);
  my $Pm = $f*(1-$e) + (1-$f)*($e/3);
  my $Pr0 = 1-$e;
  my $Pm0 = $e/3;
  my $lmut = ($Pr**$r)*($Pm**$v);
  my $lref = ($Pr0**$r)*($Pm0**$v);
  print STDERR "\t$e\t$f\t$v\t$d\t$r\t$Pr\t$Pm\t$Pr0\t$Pm0\t";
  if ($lref == 0) {
    return(1000);
  } elsif ($lmut == 0) {
    return(-1000);
  }
  my $lod = log10($lmut/$lref);
  if ($lod eq 'inf'){
    return(1000);
  } elsif ($lod eq '-inf'){
    return(-1000);
  }
  return($lod);
}

sub calNormalLOD {
  my $e = shift;
  my $f = shift;
  my $v = shift;
  my $d = shift;
  my $r = $d-$v;
  my $fg = 0.5;
  my $Pr0 = 1-$e;
  my $Pm0 = $e/3;
  my $Prg = $fg*($e/3) + (1-$fg)*(1-$e);
  my $Pmg = $fg*(1-$e) + (1-$fg)*($e/3);
  my $lref = ($Pr0**$r)*($Pm0**$v);
  my $lger = ($Prg**$r)*($Pmg**$v);
  print STDERR "\t$e\t$f\t$v\t$d\t$r\t$Pr0\t$Pm0\t$Prg\t$Pmg\t";
  if ($lref == 0) {
    return(-1000);
  } elsif ($lger == 0) {
    return(1000);
  }
  my $lod = log10($lref/$lger);
  if ($lod eq 'inf'){
    return(1000);
  } elsif ($lod eq '-inf'){
    return(-1000);
  }
  return($lod);
}

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}
