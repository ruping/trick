use strict;

my $file = shift;
my $task = shift;
my $minMark = shift;

if ($task eq '') {
  $task = "ntnb";
}
if ($minMark eq ''){
  $minMark = 20;
}

open IN, "$file";

my %colindexes;

if ($task eq 'call') {
  print "#chr\tstart\tend\tlogcopynumberratio\tallelicratio\tLOHcall\tntnb\n";
}
while ( <IN> ) {

  chomp;
  my @cols = split /\t/;
  if (/^chrom/) {
    for (my $i = 0; $i <= $#cols; $i++) {
      $colindexes{$cols[$i]} = $i;
    }
    next;
  }

  my $chrom = $cols[$colindexes{'chrom'}];
  my $start = $cols[$colindexes{'loc.start'}];
  my $end = $cols[$colindexes{'loc.end'}];
  my $numMark = $cols[$colindexes{'num.mark'}];
  my $segMean = $cols[$colindexes{'seg.mean'}];
  my $nt = $cols[$colindexes{'copynumber'}];
  my $allelicratio = $cols[$colindexes{'allelicratio'}];
  my $LOHcall = $cols[$colindexes{'LOHcall'}];
  my $cellularprevalence = $cols[$colindexes{'cellularprevalence'}];
  my $ploidy = $cols[$colindexes{'ploidy'}];
  my $normalproportion = $cols[$colindexes{'normalproportion'}];
  my $logcopynumberratio = $cols[$colindexes{'logcopynumberratio'}];

  next if ($cellularprevalence eq 'NA');

  #print STDERR "$cellularprevalence\n";
  my $nb;
  if ( exists($colindexes{'minor_cn'}) ) {
    $nb = $cols[$colindexes{'minor_cn'}];
  } else {
    $nb = round($nt*(1-(($allelicratio - $normalproportion*0.5)/(1-$normalproportion) - (1-$cellularprevalence)*0.5)/$cellularprevalence));
    $cellularprevalence = 0.99 if $cellularprevalence == 1;
  }

  next if $numMark < $minMark;
  if ($chrom =~ /^chr/) {
    $chrom =~ s/^chr/hs/;
  } else {
    $chrom = 'hs'.$chrom;
  }

  if ($task eq 'ntnb') {
    print "$chrom $start $end $nt\n";
    print STDERR "$chrom $start $end $nb\n";
  } elsif ($task eq 'raw') {
    print "$chrom $start $end $logcopynumberratio\n";
    print STDERR "$chrom $start $end $allelicratio\n";
  } elsif ($task eq 'call') {
    $chrom =~ s/^hs//;
    print "$chrom\t$start\t$end\t$logcopynumberratio\t$allelicratio\t$LOHcall\t$nt$nb\n";
  }
}

close IN;

exit 0;

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}
