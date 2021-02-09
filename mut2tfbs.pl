use strict;
use Data::Dumper;
use File::Glob ':glob';

my $mut = shift;
my $trap = shift;
my $fasta = shift;
my $pwm = shift;
my $thresdir = shift;

open PWM, "$pwm";
my $tf;
my %PWM;
my $nucl = 1;
my %nucl;
$nucl{1} = 'A';
$nucl{2} = 'C';
$nucl{3} = 'G';
$nucl{4} = 'T';
while ( <PWM> ) {
  chomp;
  if ($_ =~ /^\>\s+(\S+)/){
    $tf = $1;
  } elsif ( $_ =~ /^[\-\d]/ ) {
    my @cols = split /\s+/;
    @{$PWM{$tf}{$nucl{$nucl}}} = @cols;
    $nucl = ($nucl == 4)? 1:$nucl+1;
  } else {
    next;
  }
}
close PWM;

#print Dumper(\%PWM);


open FA, "$fasta";
my $coor;
my %fasta;
while ( <FA> ){
  chomp;
  if ($_ =~ /^\>(\S+)$/){
    $coor = $1;
  } else {
    s/\n//g;
    s/\s//g;
    $fasta{$coor} .= $_;
  }
}
close FA;

#print Dumper(\%fasta);


open TRAP, "$trap";
my %badtf;
$badtf{'ZN589_f1'}='';
$badtf{'HXD9_f1'}='';
$badtf{'DLX2_f1'}='';
$badtf{'UBIP1_f1'}='';
my %trap;
my %threshold;
my %plinear;
my %thloaded;
while ( <TRAP> ) {
  chomp;
  my ($coor, @cols) = split /\t/;
  $coor =~ s/^chr//;
  foreach my $col (@cols) {
    $col =~ /^(\S+?)\=([\d\.]+?)$/;
    my $tf = $1;
    my $score = $2;
    next if exists($badtf{$tf});

    unless ( exists($thloaded{$tf}) ) {      #load threshold
      my $threshold = "$thresdir/$tf\_thresholds.txt"; #load threshold
      my $plinear = "$thresdir/$tf\_plinear";
      if (-e $threshold) {
        open TH, "$threshold";
        while ( <TH> ) {
          chomp;
          if (/^thres/) {
            next;
          } else {
            my ($thres, $p, @others) = split /\t/;
            $threshold{$tf}{$p} = $thres;
          }
        }
        close TH;
        $thloaded{$tf} = '';

        open PL, "$plinear";
        my $pline = 1;
        while ( <PL> ){
          chomp;
          s/[\s\n]+//g;
          if ($pline == 1){
            $plinear{$tf}{'intercept'} = $_;
          } elsif ($pline == 2){
            $plinear{$tf}{'slope'} = $_;
          }
          $pline += 1;
        }
        close PL;

      } else {
        next;
      }
    }
    $trap{$coor}{$tf} = $score;
  }
}
close TRAP;

#print Dumper(\%trap);
#print "$threshold{'IKZF1_f1'}{'1.0e-05'}\n";
#print Dumper(\%threshold);
print STDERR Dumper(\%plinear);
print STDERR "tfs loaded\n";


open MUT, "$mut";
while ( <MUT> ) {
  chomp;
  if (/^[\#]?chr\t/) {
    print "$_\tmotifBreak\tmotifBreakDir\tneglog10(motifPvalue)\n";
    next;
  }
  my @cols = split /\t/;
  my $chr = $cols[0];
  my $pos = $cols[1];
  my $ref = $cols[3];
  my $alt = $cols[4];
  my $chrdhs = $cols[125];
  $chrdhs =~ s/^chr//;
  my $startdhs = $cols[126];
  my $enddhs = $cols[127];
  my $coor = $chrdhs.':'.$startdhs.'-'.$enddhs;
  my $seq = $fasta{$coor};

  my $result = 'NA';
  foreach my $tf (keys %{$trap{$coor}}) {
    next if (!exists ($PWM{$tf}));
    #motif length
    my $l = scalar(@{$PWM{$tf}{'A'}});
    #determine the scan space
    my $startscan = (($pos - $l + 1) >= $startdhs)? (($pos - $l + 1)-$startdhs+1):1;
    my $endscan   = (($pos + $l - 1) <= $enddhs)? (($pos + $l - 1)-$startdhs+1):($enddhs-$startdhs+1);
    my $scanseq = substr($seq, $startscan-1, ($endscan-$startscan+1));
    #scan the sequence to find motif
    my @scanseq = split("",$scanseq);
    my $oldSum;
    my $max = 0;
    my $maxind = 0;
    my $maxmotif = '';
    for (0..(($endscan-$l+1)-$startscan)) { #loop over all the possible motifs
      my $index = $_;
      my $sum = 0;
      my $cmotif = '';
      for (my $i = $index; $i <= ($index + $l - 1); $i++) {
        my $base = $scanseq[$i];
        $sum += ${$PWM{$tf}{$base}}[$i-$index];
        $cmotif .= $base;
      }
      if ($sum > $max) {
        $max = $sum;
        $maxind = $index;
        $maxmotif = $cmotif;
      }
    } #loop over all possible motifs

    my $thres = $threshold{$tf}{'1.0e-05'};
    if ($thres eq ''){$thres = $threshold{$tf}{'5.0e-05'}};
    if ($thres eq ''){$thres = 8};
    if ($max >= $thres) {
      $result = '';
      $result .= $tf.'=';
      my $motifStart = ($startdhs + $startscan - 1) + $maxind;
      my $mutind = $pos - $motifStart;
      my $type = (${$PWM{$tf}{$ref}}[$mutind] < ${$PWM{$tf}{$alt}}[$mutind])? '+':'-';
      $result .= $motifStart.','.$max.','.$mutind.','.$maxmotif."\t".$type;
      my $pvalue = $plinear{$tf}{'intercept'} + $plinear{$tf}{'slope'}*$max;
      $result .= "\t".$pvalue;
    }
  } #for each tf

  print "$_\t$result\n";
}
close MUT;
