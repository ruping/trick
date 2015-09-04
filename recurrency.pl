use strict;
use Data::Dumper;

my $file = shift;
my $maf = shift;
my $somaticInfo = shift; #if for somatic judgement
my $bloodCall = shift;   #whether the blood is also single called


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
  #print STDERR Dumper (\%somatic);
  #print STDERR Dumper (\%germline);
}


my @rectum = qw(AC57maf AC439maf AC440maf AC441maf AC443maf AC447maf AC525maf AC526maf AC527maf AC528maf AC529maf AC530maf AC531maf AC532maf AC533maf AC546maf AC548maf AC580maf AC637maf AC653maf AC668maf);
my @ileum = qw(AC444maf AC445maf AC446maf AC516maf AC517maf AC518maf AC519maf);
my @primary = qw(AC532maf AC533maf AC546maf AC580maf AC668maf);
my @blood = qw(AC1maf AC547maf AC581maf AC669maf);
my @all = ();
push(@all, @rectum);
push(@all, @ileum);
push(@all, "AC442maf");


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
    if ($maf eq '') {
      print "$_\tfounds\tfounds.rectum\tfounds.ileum\tfounds.primary\n";
    } elsif ($maf == 1) {
      print "$_\tmaf\n";
    } elsif ($maf =~ /trace/) {
      print "$_\ttrace\n";
    } elsif ($maf eq 'founds') {
      print "$_\tfounds\n";
    } elsif ($maf =~ /filter/) {
      print "$_\tfilter\n";
    } elsif ($maf =~ /somatic/) {
      print "$_\tsomatic\tgermline\n";
    } elsif ($maf eq 'depth'){
      print "$_\tdepthav\n";
    } elsif ($maf eq 'samfounds'){
      print "$_\tfounds\trsam\n";
    }
  } else {
    my @cols = split /\t/;
    if ($maf eq ''){
      my $founds = 0;
      my $foundsRectum = 0;
      my $foundsIleum = 0;
      my $foundsPrimary = 0;
      foreach my $rec ( @rectum ) {
        if ($cols[$colnames{$rec}] >= 0.1) {
          my $vard = sprintf("%.1f", $cols[$colnames{$rec}]*$cols[$colnames{$rec}+1]);
          if ($vard >= 2) {
            $foundsRectum++;
            $founds++;
          }
        }
      }

      foreach my $ile ( @ileum ) {
        if ($cols[$colnames{$ile}] >= 0.1) {
          my $vard = sprintf("%.1f", $cols[$colnames{$ile}]*$cols[$colnames{$ile}+1]);
          if ($vard >= 2) {
            $foundsIleum++;
            $founds++;
          }
        }
      }

      foreach my $pri (@primary) {
        if ($cols[$colnames{$pri}] >= 0.1) {
          my $vard = sprintf("%.1f", $cols[$colnames{$pri}]*$cols[$colnames{$pri}+1]);
          if ($vard >= 2) {
            $foundsPrimary++;
          }
        }
      }

      print "$_\t$founds\t$foundsRectum\t$foundsIleum\t$foundsPrimary\n";
    } elsif ($maf == 1) {    #get real maf
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
    } elsif ($maf =~ /trace/) {  #trace sample
      my $trace = '';
      foreach my $sample (@all) {
        my $samp = $sample;
        if ($samp ne 'AC3maf') {
          $samp =~ s/maf$//;
          if ($cols[$colnames{$samp}] > 0){
              $trace .= "$samp,"
          }
        } else {  #AC3maf
          if ($cols[$colnames{$samp}] >= 0.1) {
            my $vard = sprintf("%.1f", $cols[$colnames{$samp}]*$cols[$colnames{$samp}+1]);
            if ($vard >= 2) {
              $samp =~ s/maf$//;
              $trace .= "$samp,";
            }
          }
        }
      }
      $trace =~ s/,$//;
      if ($trace ne '') {
        my @trace = split(/\,/, $trace);
        my $ftrace = '';
        if (scalar(@trace) > 1) {
          if ($trace[0] eq 'AC3') {
            $ftrace = $trace[1];
          } else {
            $ftrace = $trace[0];
          }
        } else {
          $ftrace = $trace[0];
        }
        print "$_\t$ftrace\n" if ($maf eq 'trace');
        print "$_\t$trace\n" if ($maf eq 'traceall');
      }
    } elsif ($maf =~ /founds/) {    #trace all samples and check whether it is originally called
      my $founds = 0;
      my $samfounds = 0;
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colindex{$i} =~ /^(.+?)maf$/) {
          my $samp = $1;
          my $maf = $cols[$i];
          my $endsratio = 0;
          my $cmean = 0;
          my $cmedian = 0;
          my $strandRatio = 0;

          if ($cols[$i] =~ /\|/) { #split the var surrounding information
            my @infos = split(/\|/, $cols[$i]);
            $maf = $infos[0];
            $endsratio = $infos[1];
            ($cmean, $cmedian) = split(',', $infos[2]);
            $strandRatio = $infos[3];
          }

          my $depth = $cols[$i+1];
          my $vard = sprintf("%.1f", $maf*$depth);

          if (($endsratio <= 0.9 or ((1-$endsratio)*$vard >= 2)) and ($strandRatio > 0 and $strandRatio < 1) and (($cmean+$cmedian) < 5.5 or $cmedian <= 2)) {  #it looks good
            if ($maf >= 0.03 and $vard >= 3) {
              if ($somaticInfo ne '') {        #count only tumor
                if ( exists($somatic{$samp}) ) {
                  $founds++;
                  if ($colindex{$i-1} =~ /^$samp$/){  #calling information for the same sample
                    if ( $cols[$i-1] =~ /\|/ ){  #it is called originally
                      $samfounds++;
                    }
                  }
                }
              } else {
                $founds++;                     #count all samples
                if ($colindex{$i-1} =~ /^$samp$/){  #calling information for the same sample
                  if ( $cols[$i-1] =~ /\|/ ) { #it is called originally
                    $samfounds++;
                  }
                }
              }
            } #maf and vard requirements
          }  # it looks good
        } #maf
      } #each column
      print "$_\t$founds";
      if ($maf =~ /sam/){  #original founds
        my $rsam = ($founds > 0)? sprintf("%.2f", $samfounds/$founds) : 0;
        print "\t$rsam";
      }
      print "\n";
    } elsif ($maf eq 'depth') {   #find av depth
      my $dep = 0;
      my $Ndep = 0;
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colindex{$i} =~ /maf$/) {
          if ($cols[$i] >= 0.1){
            my $vard = sprintf("%.1f", $cols[$i]*$cols[$i+1]);
            if ($vard >= 2) {
              $Ndep++;
              $dep += $cols[$i+1];
            }
          }
        } #maf
      } #each column
      my $depav = ($Ndep > 0)? sprintf("%.1f", $dep/$Ndep):0;
      print "$_\t$depav\n";
    } elsif ($maf =~ /filter/) {     #filter (based on Jie's Results)
      my @detectedSample;
      my %detectedSample;
      my %mafs;
      my $chr;
      my $pos;
      my $endsratio = 0;
      my $strandRatio = 0;
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
            my @traceSomatic = split(',', $traceSomatic);
            push(@detectedSample, @traceSomatic);
          }
          if ($traceGermline =~ /\,/) {
            my @traceGermline = split(',', $traceGermline);
            push(@detectedSample, @traceGermline);
          }
          foreach my $detectedSamp (@detectedSample) {
            $detectedSample{$detectedSamp} = '';
          }
          #print STDERR "$traceSomatic\t$traceGermline\n";
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
          $mmaf = ($infos[0] > $mmaf)? $infos[0]:$mmaf;
          $strandRatio = ($infos[0] > $mmaf)? $infos[3]:$strandRatio;
        }
      }

      my $status;
      print STDERR "$chr\t$pos\t$rep$sc\t$detectedSample[0]\t$mmaf\t$endsratio\t$cmean\t$cmedian\t$cmeanav\t$cmedianav\n";
      if ($rep == 1 and $sc == 1) {
        $status = ($endsratio < 0.9 and ($strandRatio > 0 and $strandRatio < 1) and (($cmean+$cmedian) < 4.5 or $cmedian < 2) and ($cmeanav + $cmedianav) < 5.2)? 'PASS':'FOUT';   #conservative for rep and sc
      } elsif ($rep == 1 or $sc == 1) {
        $status = ($endsratio < 0.9 and ($strandRatio > 0 and $strandRatio < 1) and (($cmean+$cmedian) < 5 or $cmedian <= 2) and ($cmeanav + $cmedianav) < 5.2)? 'PASS':'FOUT';
      } else {
        $status = ($endsratio < 0.9 and ($strandRatio > 0 and $strandRatio < 1) and (($cmean+$cmedian) < 5.5 or $cmedian <= 2) and ($cmeanav + $cmedianav) < 5.5)? 'PASS':'FOUT';
      }
      print "$_\t$status\n" if ($status eq 'PASS');
    } elsif ($maf =~ /somatic/) {  #find somatic ones
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
          my $cmean = 0;
          my $cmedian = 0;

          if ($cols[$i] =~ /\|/) {  #split the var surrounding information
            my @infos = split(/\|/, $cols[$i]);
            $maf = $infos[0];
            $endsratio = $infos[1];
            ($cmean, $cmedian) = split(',', $infos[2]);
            $strandRatio = $infos[3];
          }

          my $depth = $cols[$i+1];
          my $vard = sprintf("%.1f", $maf*$depth);

          if (exists $somatic{$samp}) {     #for tumor samples require some additional thing
            if (($endsratio <= 0.9 or ((1-$endsratio)*$vard >= 2)) and ($strandRatio > 0 and $strandRatio < 1) and (($cmean+$cmedian) < 5.2 or $cmedian <= 2)) { #true event
              $maf = $maf;
            } else {
              $maf = 0;         #not reliable somatic
            }
          }

          #print STDERR "$samp\t$maf\t$endsratio\t$cmean\t$cmedian\n";

          if (exists $somatic{$samp}) {       #it is tumor
             $tumor{$samp} = $maf if ($vard >= 3 and $maf >= 0.03);
          } elsif (exists $germline{$samp}) { #it is blood
            foreach my $ct (@{$germline{$samp}}) {
              if ($bloodCall eq 'yes' and $cols[$i-1] =~ /\|/) {       #it is originally called
                $bloodcalled{$ct} = '';
              }
              if ( $maf == 0 and $depth >= 8 ) {
                $nonblood{$ct} = '';
              } elsif ( $maf == 0 and $depth < 8 ) {
                $unknown{$ct} = '';
              } else {                          #maf != 0
                $blood{$ct} = $maf;
              }
            }
          } #it is blood
        } #maf
      } #each column

      my $soma = 'NA';
      my $germ = 'NA';
      foreach my $tumorSamp (keys %tumor) {
        my $stype = 'NA';
        if (exists $nonblood{$tumorSamp}) {
          $stype = 'good';
          $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';
        } elsif (exists($blood{$tumorSamp})) {
          if ($blood{$tumorSamp} < 0.03 and $tumor{$tumorSamp}/$blood{$tumorSamp} >= 4) {
            $stype = 'doubt';
            $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';
          } elsif (exists($unknown{$tumorSamp})) {
            $stype = 'undef';
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
