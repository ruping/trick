use strict;
use lib '/ifs/home/c2b2/ac_lab/rs3412/cpan/lib/perl5/site_perl/5.10.0';
use Statistics::Basic qw(:all nofill);

my $fusion = shift;
my $sample = shift;
my $task = shift;

my %need;
foreach my $need ( split(',',$sample) ) {
  $need{$need} = '';
}

my %fusions;

open IN, "$fusion";
#my $id = 1;
while ( <IN> ){
  chomp;
  my ($gene1, $gene2, $bp1, $bp2, $dir, $rep, $sc, $type, $strand, $cov1, $cov2, $cov3, $cov4, $cov5, $transcripts) = split /\t/;

  my $sample;
  my $g1;
  my $g2;
  if ($gene1 =~ /(AC\d+)\.fusion\.report\:(\w+?)\((.+?)\)$/) {
    $sample = $1;
    next if (!exists($need{$sample}));
    $g1 = $3;
  } elsif ($gene1 =~ /(AC\d+)\.fusion\.report\:IGR$/) {
    $sample = $1;
    next if (!exists($need{$sample}));
    $g1 = 'IGR';
  } else {
    print STDERR "wierd gene name: $gene1\n";
    exit 22;
  }


  if ($gene2 =~ /^(\w+?)\((.+?)\)$/) {
    $g2 = $2;
  } elsif ($gene2 eq 'IGR') {
    $g2 = 'IGR';
  } else {
    print STDERR "wierd gene name: $gene2\n";
    exit 22;
  }


  if (($g2 ne 'IGR' and $g1 =~ /$g2/) or ($g1 ne 'IGR' and $g2 =~ /$g1/)) {
     next;
  }
  if (($g1 eq 'TMPRSS3' and $g2 eq 'DSCAM') or ($g2 eq 'TMPRSS3' and $g1 eq 'DSCAM')) {
     next;
  }
  if (($g1 =~ /^IGK[JV]/ or $g2 =~ /^IGK[JV]/) and $sc eq 'CC') {
     next;
  }
  if (($g1 eq 'IGR' and $g2 =~ /^RP\d+\-\d+/) or ($g2 eq 'IGR' and $g1 =~ /^RP\d+\-\d+/)){
     next;
  }
  if (($g1 =~ /^$g2/ or $g2 =~ /^$g1/) and $sc eq 'CC') {
     next;
  }
  if ($rep =~ /R/ and $sc =~ /C/) {
     next;
  }
  if ($cov5 <= 1.5){
     next;
  }
  if ($cov5 < 2 and $sc eq 'CC'){
     next;
  }
  if ($cov4/$cov5 >= 10 and $cov5 <= 3) {
     next;
  }
  if ($cov4/$cov5 >= 20 and $cov5 < 10) {
     next;
  }
  if ($sc eq 'CC') {
     next;
  }
  if ($cov3 < 2){
     next;
  }


  my $options = "thickness=$cov5";


  $bp1 =~ /^(.+?)\:(\d+)$/;
  my $chr1 = $1;
  my $pos1 = $2;
  $chr1 =~ s/^chr//;
  $chr1 = 'hs'.$chr1;
  $bp2 =~ /^(.+?)\:(\d+)$/;
  my $chr2 = $1;
  my $pos2 = $2;
  $chr2 =~ s/^chr//;
  $chr2 = 'hs'.$chr2;

  if ($chr1 eq $chr2 and abs($pos1 - $pos2) < 500000){
      next;
  }

  my $key = join("\t", sort {$a cmp $b} ($g1, $g2));
  $fusions{$key}{$sample}{'support'} = $cov5;
  $fusions{$key}{$sample}{'expression'} = $cov4;
  $fusions{$key}{$sample}{'coor'} = "$chr1 $pos1 $pos1 $chr2 $pos2 $pos2";

  if ($task ne 'merge'){

    if ($task eq 'original'){
      print "$_\n";
      next;
    }

    print "$chr1 $pos1 $pos1 ";
    print "$chr2 $pos2 $pos2 $options\n";

  }
}
close IN;


if ( $task eq 'merge' ) {
  foreach my $fusion (keys %fusions){
    my $samps;
    my @supports;
    my @expressions;
    my $coor;
    foreach my $samp (keys %{$fusions{$fusion}}){
      $samps .= $samp."\t";
      push(@supports, $fusions{$fusion}{$samp}{'support'});
      push(@expressions, $fusions{$fusion}{$samp}{'expression'});
      if ($coor eq ''){
        $coor = $fusions{$fusion}{$samp}{'coor'};
      }
    }
    my $recur = scalar(keys %{$fusions{$fusion}});
    my $support = sprintf("%.1f", mean(@supports));
    my $expression = sprintf("%.1f", mean(@expressions));
    print STDERR "$fusion\t$recur\t$samps\t$support\t$expression\n";
    my $rsup = round($support);
    #$rsup .= 'p,';
    #my $options = "thickness=$rsup";
    my @coors = split (/\s/, $coor);
    $coors[1] -= ($coors[1] > ($rsup*800000))? $rsup*800000 : ($coors[1]-1);
    $coors[2] += $rsup*800000;
    $coors[4] -= ($coors[4] > ($rsup*800000))? $rsup*800000 : ($coors[4]-1);
    $coors[5] += $rsup*800000;
    $coor = join(' ', @coors);
    my $options;
    if ($recur > 6) {
      $options .= "color=dturquoisetrans";
    } elsif ($recur >= 2) {
      $options .= "color=turquoisetrans";
    } elsif ($recur == 1) {
      $options .= "color=vwblacktrans";
    } else {
      $options .= "color=vvwblacktrans";
    }

    print "$coor $options\n";
  }

} #merge for circos


sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}
