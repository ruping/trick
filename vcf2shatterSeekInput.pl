use Data::Dumper;
use strict;

my $vcf = shift;
my $sample = shift;
my $split = shift;
if ($split eq '') {
  $split = "no";
}
#print STDERR "task: split = $split\n";

open IN, "$vcf";
my %colnames;
my %colindex;
if ($split eq 'no') {
  printf("%s\n", join("\t", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "sv_id", "pe_support", "strand1", "strand2", "svclass", "svmethod"));
} elsif ($split eq 'vcf') {
} else {
  printf("%s\n", join("\t", "chrom", "start", "end", "sv_id", "pe_support", "strand", "svclass", "svmethod"));
}

while ( <IN> ) {

  chomp;
  if (/^\#\#/){
    print "$_\n" if ($split eq 'vcf');
    next;
  }

  if (/^#CHROM/) {  #it is header
    $_ =~ s/^#//;
    my @cols = split /\t/;
    for(my $i = 0; $i <= $#cols; $i++) {
      $colnames{$cols[$i]} = $i;
      $colindex{$i} = $cols[$i];
    }
    if ($split eq 'vcf'){
      $cols[0] = '#'.$cols[0];
      printf("%s\n", join("\t", @cols[0..8], $sample));
    }

  } else {        #variants
    my @cols = split /\t/;
    my @formats = split(':', $cols[$colnames{'FORMAT'}]);
    my %formindex;
    for (my $f = 0; $f <= $#formats; $f++) {
      $formindex{$formats[$f]} = $f;
    }
    my @sdata = split(':', $cols[$colnames{$sample}]);

    if ($sdata[$formindex{'GT'}] !~ /1/) {        #genotyping do not contain 1
      next;
    }

    my $info = $cols[$colnames{'INFO'}];
    my $sv_id = $cols[$colnames{'ID'}];

    #processing
    (my $chrom1 = $cols[$colnames{'CHROM'}]) =~ s/^chr//;
    my $pos1 = $cols[$colnames{'POS'}];
    my $start1;
    my $end1;
    my $chrom2;
    my $pos2;
    my $start2;
    my $end2;
    my $pe_support;
    my $svclass;
    my $strand1;
    my $strand2;
    my $svmethod;
    my $invOrt;

    if ($info =~ /CHR2\=([a-z0-9XYMT]+);/) {
      ($chrom2 = $1) =~ s/^chr//;
    }
    if ($info =~ /END\=(\d+);/) {
      $pos2 = $1;
    }
    if ($info =~ /CIPOS\=(-\d+)\,(\d+);/) {
      $start1 = $pos1 + $1;
      $end1 = $pos1 + $2;
    }
    if ($info =~ /CIEND\=(-\d+)\,(\d+);/) {
      $start2 = $pos2 + $1;
      $end2 = $pos2 + $2;
    }
    if ($info =~ /PE\=(\d+);/) {
      $pe_support = $1;
    }
    if ($info =~ /SVTYPE=([A-Z]+);/) {
      $svclass = $1;
      $svclass =~ s/BND/TRA/;
    }
    if ($info =~ /SVMETHOD=(.+?)\;/) {
      $svmethod = $1;
    }
    my $ct;
    if ($info =~ /CT=(([35])to([35]))\;/) {
      $ct = $1;
      my $ort1 = $2;
      my $ort2 = $3;
      ($invOrt = $ct) =~ s/to/2/;
      $invOrt =~ tr/35/th/;
      $svclass = ($svclass eq "INV")? $invOrt.$svclass : $svclass;

      $strand1 = ($ort1 == 3)? '+' : '-';
      $strand2 = ($ort2 == 3)? '+' : '-';
    }
    next if ($chrom1 eq 'Y' or $chrom2 eq 'Y');  #no Y chromosome
    if ($split eq 'no') {
      printf("%s\n", join("\t", $chrom1, $start1, $end1, $chrom2, $start2, $end2, $sv_id, $pe_support, $strand1, $strand2, $svclass, $svmethod));
    } elsif ($split eq 'vcf') {
      #print STDERR Dumper(\@cols);
      printf("%s\n", join("\t", @cols[0..8], $cols[$colnames{$sample}]));
    } else {
      printf("%s\n", join("\t", $chrom1, $start1, $end1, $sv_id, $pe_support, $strand1, $svclass, $svmethod));
      printf("%s\n", join("\t", $chrom2, $start2, $end2, $sv_id, $pe_support, $strand2, $svclass, $svmethod));
    }
  }

}
close IN;
