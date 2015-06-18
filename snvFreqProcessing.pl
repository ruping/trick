use strict;

my $file = shift;
my $commonAsSomatic = shift;
my %commonAsSomatic;
my %somaticTotal;

open IN, "$file";
my %colindex;
while ( <IN> ) {
  chomp;
  my @cols = split /\t/;
  if ($_ =~ /^[\#]?chr\t/){  #header rows
    for (my $i = 0; $i <= $#cols; $i++) {
      $colindex{$cols[$i]} = $i;
    }
    print "$_\tfreq\n" if ($commonAsSomatic eq '');
  } else {  #data rows
    my $id = $cols[$colindex{'id'}];
    my $esp = $cols[$colindex{'ESP6500SI.snv.vcf.INFO.MAF'}];
    my $dbsnp = $cols[$colindex{'dbsnp142.snv.vcf.INFO.CAF'}];
    my $somatic = $cols[$colindex{'somatic'}];
    my $founds = $cols[$colindex{'founds'}];
    my $espTac = $cols[$colindex{'ESP6500SI.snv.vcf.INFO.TAC'}];
    my $freq = 'NA';

    my $keep = 0;
    # somatic Ones
    if ($somatic ne 'NA'){
      $keep = 1;
    }

    #unkonwn ones
    if ($id eq '.' and $dbsnp eq 'NA' && $esp eq 'NA') {
      $keep = 1;
    }

    # known DBsnp rare ones
    my $freqDBsnp = 22;
    if ($dbsnp =~ /^\[(.+?)\]$/) {
      my @dbf = split(",", $1);
      $freqDBsnp = $dbf[$#dbf];
      if ($freqDBsnp < 0.005) {
        $keep = 1;
      }
    }

    # known ESP rare ones
    my $freqESP = 22;
    if ($freqDBsnp == 22 and $esp ne 'NA' and $espTac ne '') {
      my @espf = split(",", $esp);
      my @espTac = split(",", $espTac);  #in the format of alt, ref
      if ($espTac[0] < $espTac[1]) {       #alt should be minor allele
        $freqESP = $espf[$#espf] * 0.01;
      }
      if ($freqESP < 0.005) {
        $keep = 1;
      }
    }

    my $freqESP5400 = 22;
    if ($freqDBsnp == 22 and $freqESP == 22 and $id =~ /ESP5400\=(.+)/) {
      $freqESP5400 = $1;
      if ($freqESP5400 < 0.005) {
        $keep = 1;
      }
    }


    if ($freqDBsnp != 22) {
      $freq = $freqDBsnp;
    } elsif ($freqESP != 22) {
      $freq = $freqESP;
    } elsif ($freqESP5400 != 22){
      $freq = $freqESP5400;
    }

    if ($commonAsSomatic eq 'cs') { #commonAsSomatic
      next if ($somatic eq 'NA');
      my @somatic = split(",", $somatic);
      foreach my $ssamp (@somatic) {
        $ssamp =~ /^(.+?)\[/;
        my $ss = $1;
        $somaticTotal{$ss}++;
      }
      next if ($freq eq 'NA' or $freq < 0.005);
      #now only the somatic ones with MAF > 0.005 are retained
      foreach my $ssamp (@somatic) {
        $ssamp =~ /^(.+?)\[/;
        my $ss = $1;
        $commonAsSomatic{$freq}{$ss}++;
      }
      next;
    }

    if ($keep == 1){
      print "$_\t$freq\n";
    }
  }
}
close IN;


if ($commonAsSomatic eq 'cs'){
  print "freq";
  foreach my $ssamp (sort {$a =~ /(\d+)/; my $da = $1; $b =~ /(\d+)/; my $db = $1; $da <=> $db} keys %somaticTotal){
    print "\t$ssamp";
  }
  print "\n";
  foreach my $freq (sort {$a <=> $b} keys %commonAsSomatic) {
    print "$freq";
    foreach my $ssamp (sort {$a =~ /(\d+)/; my $da = $1; $b =~ /(\d+)/; my $db = $1; $da <=> $db} keys %somaticTotal) {
      my $frac;
      if ($commonAsSomatic{$freq}{$ssamp} eq '' or $somaticTotal{$ssamp} == 0){
        $frac = 0;
      } else {
        $frac = sprintf("%.3f", $commonAsSomatic{$freq}{$ssamp}/$somaticTotal{$ssamp});
      }
      print "\t$frac";
    }
    print "\n";
  }
}
