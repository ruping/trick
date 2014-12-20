use strict;
use File::Glob ':glob';

my $list = shift;   #filename of all vcfs
my $type = shift;
my $recheck = shift;
my $task = shift;

my @list;
open IN, "$list";
while ( <IN> ) {
  chomp;
  next if /^#/;
  push(@list, $_);
}
close IN;

open DR, "/ifs/scratch/c2b2/ac_lab/rs3412/no1/net/dna2rna.mapping";
my %rna2dna;
while ( <DR> ) {
  #dna     rna     availability
  next if /^dna\trna/;
  my ($dna, $rna, $availability) = split /\t/;
  if ($availability eq 'N') {
    next;
  }
  $rna2dna{$rna} = $dna;
}
close DR;

my %somatic;
my %samples;
foreach my $file (@list) {
  my $name;
  if ($task ne 'tcga' and $file =~ /(AC\d+)[^a-zA-Z0-9]/) {
    $name = $1;
  }
  if ($task eq 'tcga') {
    $file =~ /\/((TCGA\-[^\-]+\-[^\-]+)\-[^\-]+\-[^\-]+\-[^\-]+\-\d+)\./;
    $name = $1;
  }
  print STDERR "$name\t$file\n";

  if ($task eq 'rnaediting'){
    $name = $rna2dna{$name};
  }

  $samples{$name} = '';

  my %recheck;
  if ($recheck ne '') {
     open RECHECK, "$recheck/$type/$name";
     while ( <RECHECK> ) {
        chomp;
        my @cols = split /\t/;
        my $coor = $cols[0].':'.$cols[1];
        if ($type eq 'snv'){
          if (($cols[4]+$cols[5]+$cols[6]+$cols[7]) == 0 and $cols[2] >= 10) {                       #only record for absent indels
             $recheck{$coor} = 1;
          }
        } elsif ($type eq 'indel') {
          if ($cols[6] == 0 and $cols[5] >= 10) {                       #only record for absent indels
              $recheck{$coor} = 1;
          }
        }
     }
     close RECHECK;
     print STDERR "$recheck/$type/$name loaded\n";
  }
  my $tmpc = scalar(keys %recheck);
  print STDERR "recheced somatic: $tmpc\n";


  my $openway = $file;
  if ($file =~ /\.gz$/){
    $openway = "gzip -dc $file |";
  }
  print STDERR "$openway\n";

  open IN, "$openway";
  while ( <IN> ) {
     chomp;
     next if /^#/;
     my ($chr, $pos, $id, $ref, $alt, $qual, $pass, $info, $format, $sample, $blood) = split /\t/;

     #if ($name eq 'GT001') {
     #   my $tmp = $sample;
     #   $sample = $blood;
     #   $blood = $tmp;
     #}

     next if ($qual < 30 and $task ne 'rnaediting');
     next if ($qual < 30 and $task eq 'rnaediting');

     my @sample = split(/\:/,$sample);

     if ($sample[0] eq '0/0') {   #skip some wierd thing
        next;
     }

     ###########################################################################decide somatic
     my $somatic = 0;
     if ($type eq 'snv') {   #for snp
        my @blood = split(/\:/,$blood);
        my @bad = split (/\,/, $blood[2]);
        if ($blood[0] eq '0/0' and $bad[1] == 0) {
           $somatic = 1;
        }
        if ($recheck ne '') {       #recheck
           my $coor = $chr.':'.$pos;
           if (exists($recheck{$coor})) {
              $somatic = 1;
           }
        }
     } elsif ($type eq 'indel') {
         my @blood = split(/\:/,$blood);
         if ($blood[0] eq '0/0') {
           $somatic = 1;
         }
         if ($recheck ne ''){        #recheck
           my $coor = $chr.':'.$pos;
           if (exists($recheck{$coor})){
              $somatic = 1;
           }
         }
     }
     #############################################################################decide somatic


     next if ($task eq 'rnaediting' and $somatic == 0);                  #for RNAediting study, only remember somatic ones


     if ($task ne 'rnaediting' and ($id ne '.' or $info =~ /dbSNP/ or $info =~ /1KG\=/ or $info =~ /ESP5400\=/)) {  #snp in population, is it a somatic one?

       my $freq = -1;
       if ($info =~ /(1KG=(.+?));/) {
           my $kid = $1;                               #re define $id to 1KG when absent
           $freq = $2;
           if ($id eq '.') {
              $id = $kid;
           }
       }

       if ($info =~ /(ESP5400=(.+?));/) {
           my $kid = $1;                               #re define $id to ESP when absent
           $freq = $2;
           if ($id eq '.') {
              $id = $kid;
           }
       }

       if ($somatic == 0) {                             #keep somatic ones even if it is marked as a common snp
         if ($freq == -1) {
            next;
         }
         elsif ($freq > 0 and ($id !~ /^1KG/ and $id !~ /^ESP5400/)) {   #still dbSNP ones
            next;
         }
       }
     } #somatic common snp


     if ($info =~ /MQ0Fraction=(.+?);/) {
       next if ($1 > 0.1);
     }
     if ($info =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/) {
       my $tailb = $4;
       if ($tailb =~ /e/) {
          next;
       } elsif ($tailb < 0.005){
          next;
       }
     }

     my $coor = $chr.':'.$pos;

     my $maf = -1;
     if ($info =~ /DP4\=(\d+)\,(\d+)\,(\d+)\,(\d+)\;/) {
        $maf = sprintf("%.3f", ($3+$4)/($1+$2+$3+$4));
     } else {
        my @sampleinfo = split(":", $sample);
        $sampleinfo[2] =~ /^(\d+)\,(\d+)$/;
        $maf = sprintf("%.3f", $2/($1+$2));
     }

     $somatic{$coor}{$name} = $maf;
     my $function;
     $info =~ /(function=.+?$)/;
     $function = $1;
     $somatic{$coor}{'function'} = $function;
     $somatic{$coor}{'info'} = join("\t", ($id,$ref,$alt));
     #if ($somatic{$coor}{'somatic'} eq '') {
     #  $somatic{$coor}{'somatic'} = 0;
     #  $somatic{$coor}{'germline'} = 0;
     #}
     if ($somatic == 1) {
        #$somatic{$coor}{'somatic'} = '' if ($somatic{$coor}{'somatic'} == 0);
        $somatic{$coor}{'somatic'} .= $name.',';
     } else {
        #$somatic{$coor}{'germline'} = '' if ($somatic{$coor}{'germline'} == 0);
        $somatic{$coor}{'germline'} .= $name.',';
     }

  }
  close IN;
  $tmpc = scalar(keys %somatic);
  print STDERR "table tmp: $tmpc\n";

}

print "#chr\tpos\tid\tref\talt";
foreach my $name (sort keys %samples) {
   print "\t$name";
}
print "\tfunction\tsomatic\tgermline\n";

foreach my $coor (sort {$a =~ /^(\w+):(\d+)$/; my $ca = $1; my $pa = $2; $b =~ /^(\w+):(\d+)$/; my $cb = $1; my $pb = $2; $ca cmp $cb || $pa <=> $pb} keys %somatic) {
  $coor =~ /^(\w+):(\d+)$/;
  my $info = $somatic{$coor}{'info'};
  print "$1\t$2\t$info";
  foreach my $name (sort keys %samples) {
    if ($somatic{$coor}{$name} ne '') {
       print "\t$somatic{$coor}{$name}";
    } else {
       print "\t0";
    }
  }
  my $function = $somatic{$coor}{'function'};
  my $somatic = ($somatic{$coor}{'somatic'} eq '')? 0 : $somatic{$coor}{'somatic'};
  my $germline = ($somatic{$coor}{'germline'} eq '')? 0 : $somatic{$coor}{'germline'};
  print "\t$function\t$somatic\t$germline\n";
}
