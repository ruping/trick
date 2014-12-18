use strict;
use File::Glob ':glob';

my $list = shift;   #filename of all vcfs
my $type = shift;


my @list;
open IN, "$list";
while ( <IN> ) {
  chomp;
  next if /^#/;
  push(@list, $_);
}
close IN;

my %somatic;
my %samples;
foreach my $file (@list) {
  my $name;
  if ($file =~ /(AC\d+)[^a-zA-Z0-9]/){
    $name = $1;
  }
  print STDERR "$name\t$file\n";

  $samples{$name} = '';

  open IN, "$file";
  while ( <IN> ) {
     chomp;
     next if /^#/;
     my ($chr, $pos, $id, $ref, $alt, $qual, $pass, $info, $format, $sample, $blood) = split /\t/;

     if ($qual ne '.'){
        next if ($qual < 30);
     }

     my @sample = split(/\:/,$sample);

     if ($sample[0] eq '0/0') {   #skip some wierd thing
        next;
     }


     if ($id ne '.' or $info =~ /dbSNP/ or $info =~ /1KG/) {  #snp in population, is it a somatic one?
       my $somatic = 0;

       my $freq = -1;
       if ($info =~ /(1KG=(.+?));/) {
           my $kid = $1;                               #re define $id when absent
           $freq = $2;
           if ($id eq '.') {
              $id = $kid;
           }
       }

       if ($type eq 'snv') {   #for snp
         my @blood = split(/\:/,$blood);
         my @bad = split (/\,/, $blood[2]);
         if ($blood[0] eq '0/0' and $bad[1] == 0) {
           $somatic = 1;
         }
       } elsif ($type eq 'indel'){
         my @blood = split(/\:/,$blood);
         if ($blood[0] eq '0/0') {
           $somatic = 1;
         }
       }

       if ($somatic == 1){
          next;
       }
     } else {
        next;    #not in population, skip
     }


     if ($info =~ /MQ0Fraction=(.+?);/) {
       next if ($1 > 0.1);
     }

     if ($info =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/) {
       my $tailb = $4;
       if ($tailb =~ /e/) {
          next;
       } elsif ($tailb < 0.005) {
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
  }
  close IN;
}

print "#chr\tpos\tid\tref\talt";
foreach my $name (sort keys %samples) {
   print "\t$name";
}
print "\tfunction\n";

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
  print "\t$function\n";
}
