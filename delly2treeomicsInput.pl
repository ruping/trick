use strict;

my $delly = shift;
my $sampleRename = shift;
my $pvaf = shift;
my $skipLC = shift;
if ($pvaf eq '') {
  $pvaf = "FALSE";
}
if ($skipLC eq '') {
  $skipLC = "TRUE";
}

my %sampleRename;
open IN, "$sampleRename";
while ( <IN> ) {
  chomp;
  my ($sn, $mapping) = split("\t", $_);
  $sampleRename{$sn} = $mapping;
}
close IN;

open IN, "$delly";
my %colnames;
my %colindex;
while ( <IN> ) {
  chomp;
  next if /^##/;
  if ($_ =~ /^#CHROM/) {
    s/^#//;
    my @cols = split("\t", $_);
    for( my $i = 0; $i <= $#cols; $i++ ) {
      $colnames{$cols[$i]} = $i;
      $colindex{$i} = $cols[$i];
    }
    my @samples;
    for (my $i = 9; $i < $#cols; $i++) {  #all samples
      #my $sample = $sampleRename{$colindex{$i}};
      my $sample = $cols[$i];
      if ($pvaf eq "TRUE") {
        $sample .= "mafc\t".$sample."refc\t".$sample."altc\t".$sample."d\t".$sample."mafa\t".$sample."ccf\t".$sample."ccfSD\t".$sample."time";
      }
      push(@samples, $sample);
    }
    my @header = qw(chr pos Change Gene);
    push(@header, @samples);
    printf("%s\n", join("\t", @header));
    if ($pvaf eq "FALSE"){
      printf STDERR ("%s\n", join("\t", @header));
    }
  } else {
    my @cols = split("\t", $_);
    my $chrom = $cols[$colnames{'CHROM'}];
    my $pos = $cols[$colnames{'POS'}];
    (my $change = $cols[$colnames{'ALT'}]) =~ s/[<>]//g;
    if ( $change =~ /chr/ ) {
      $change = 'TRA';
    }
    #$change = $cols[$colnames{'REF'}].'>A';
    #$change =~ s/N/A/;
    my $gene = $cols[$colnames{'ID'}];

    #IN: CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
    my @formats = split(':', $cols[$colnames{'FORMAT'}]);
    my %formindex;
    for (my $f = 0; $f <= $#formats; $f++) {
      $formindex{$formats[$f]} = $f;
    }

    my @altc;
    my @depth;
    my @vaf;
    for (my $i = 9; $i < $#cols; $i++) {  #all samples
      my @sampleInfo = split(":", $cols[$i]);
      my $DR = $sampleInfo[$formindex{'DR'}];     #high-quality reference pairs
      my $DV = $sampleInfo[$formindex{'DV'}];     #high-quality variant pairs
      my $RR = $sampleInfo[$formindex{'RR'}];     #high-quality reference junction reads
      my $RV = $sampleInfo[$formindex{'RV'}];     #high-quality variant junction reads
      my $altc = $DV+$RV;
      my $depth = $DV+$RV+$DR+$RR;
      my $vaf = ($depth > 0)? sprintf("%.4f", $altc/$depth) : 0;
      push(@altc, $altc);
      push(@depth, $depth);
      push(@vaf, $vaf, $depth-$altc, $altc, $depth, 0, 0, 0, 0);
    }
    if ($pvaf eq "TRUE") {
      printf("%s\n", join("\t", $chrom, $pos, $change, $gene, @vaf));
    } else {
      printf("%s\n", join("\t", $chrom, $pos, $change, $gene, @altc));
      printf STDERR ("%s\n", join("\t", $chrom, $pos, $change, $gene, @depth));
    }
  }
}
close IN;

