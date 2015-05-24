use strict;
use File::Glob ':glob';
use List::Util qw[min max];
use Data::Dumper;
use Getopt::Long;

my $file;      #filename of all rechecked files
my $type;
my $original;
my $task;      #only for tcga
my $prefix;
my $blood;

GetOptions (
           "file|f=s"       => \$file,             #filename
           "type|t=s"       => \$type,             #snv or indel
           "original|o=s"   => \$original,         #original big table
           "task|k=s"       => \$task,             #task type
           "prefix|p=s"     => \$prefix,
           "blood|b=s"      => \$blood,            #blood
           "help|h"         => sub{
                               print "usage: $0 get all minor allele frequency for samples under recheck\n\nOptions:\n\t--file\t\tthe filename of all rechecked files\n";
                               print "\t--type\t\tthe type of variants, snv or indel\n";
                               print "\t--original\tthe original mutation big table\n";
                               print "\t--prefix\tthe prefix of samples' names\n";
                               print "\t--blood\tthe sample names of blood samples\n";
                               print "\t--task\t\tthe task, such as tcga\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );



my %blood;
foreach my $bl (split(/\,/, $blood)){
  $blood{$bl} = '';
}

my @list;
open IN, "$file";
while ( <IN> ) {
  chomp;
  next if /^#/;
  push(@list, $_) unless ($_ eq '');
}
close IN;


####################################################get chr pos
my %chrJumper;
$chrJumper{'original'} = getchrpos($original);

my %samples;
foreach my $file (@list) {
  my $name;
  if ($prefix ne '' and $file =~ /($prefix\d+)$/) {
    $name = $1;
  } elsif ($task eq 'tcga' and $file =~ /\/((TCGA\-([^\-]+\-[^\-]+))\-[^\-]+\-[^\-]+\-[^\-]+\-\d+)$/) {
    $name = $1;
  }
  if ($name ne '') {
    $samples{$name} = '';
    $chrJumper{$name} = getchrpos($file);
  }
}

print STDERR Dumper(\%chrJumper);
####################################################get chr pos


print "#chr\tpos\tid\tref\talt";
foreach my $name (sort {$a =~ /$prefix(\d+)/; my $ia = $1; $b =~ /$prefix(\d+)/; my $ib = $1; $ia <=> $ib} keys %samples) {
   print "\t$name\t$name".'d';
}
print "\tcmeanav\tcmedianav\n";

foreach my $chrc (sort keys %{$chrJumper{'original'}}) {

  my $jumperO = $chrJumper{'original'}->{$chrc};
  open OR, "$original";
  seek(OR, $jumperO, 0);
  my %OR;
  my $diji = 1;
  my $oldCoor = 'SRP';
  while ( <OR> ) {
    next if /^[@#]/;
    chomp;
    my @cols = split /\t/;
    my $chr = $cols[0];

    last if $chr ne $chrc;

    my $pos = $cols[1];
    my $id  = $cols[2];
    my $coor = $chr.':'.$pos;
    if ($coor eq $oldCoor) {
      $diji += 1;
    } else {
      $diji = 1;
    }
    my $ref = $cols[3];
    my $alt = $cols[4];
    $OR{$coor}{$diji} = $id.','.$ref.','.$alt;

    #redefine
    $oldCoor = $coor;
  }
  close OR;
  print STDERR "$original $chrc loaded\n";


  my %somatic;
  foreach my $file (@list) {
    my $name;
    if ($prefix ne '' and $file =~ /($prefix\d+)$/) {
      $name = $1;
    } elsif ($task eq 'tcga' and $file =~ /\/((TCGA\-([^\-]+\-[^\-]+))\-[^\-]+\-[^\-]+\-[^\-]+\-\d+)$/) {
      $name = $1;
    }

    print STDERR "$name\n";

    my $jumperI = $chrJumper{$name}->{$chrc};
    my $djindex = 1;
    my $prevCoor = 'SRP';
    open IN, "$file";
    seek(IN, $jumperI, 0);
    while ( <IN> ) {
      chomp;
      next if /^[#@]/;
      next if /^chr\t/;
      if ($type =~ /indel/) {       #indel
        my ($chr, $pos, $ref, $alt, $indelType, $depth, $vard, $junction, $cmean, $cmedian) = split /\t/;
        if ($cmean =~ /e/){
          $cmean = 0;
        }
        last if $chr ne $chrc;

        my $coor = $chr.':'.$pos;
        if ($coor eq $prevCoor) {
          $djindex += 1;
        } else {
          $djindex = 1;
        }

        if ($vard > 0 and $depth > 0) {
          $somatic{$coor}{$djindex}{$name} = sprintf("%.3f", $vard/$depth);
        } else {
          $somatic{$coor}{$djindex}{$name} = 0;
        }
        $somatic{$coor}{$djindex}{$name} .= "\t$depth";
        if ($junction != 0) {  #there are some junction reads
          $somatic{$coor}{$djindex}{$name} .= ",$junction";
        }
        $somatic{$coor}{$djindex}{'info'} = join("\t", ($ref,$alt));
        $somatic{$coor}{$djindex}{'consecutive'} .= $cmean.','.$cmedian.','.$vard.';';

        #redefine
        $prevCoor = $coor;

      } elsif ($type =~ /snv/) {    #snv

        my ($chr, $pos, $depth, $vard, $A, $C, $G, $T, $vends, $junction, $cmean, $cmedian) = split /\t/;
        if ($cmean =~ /e/) {
          $cmean = 0;
        }
        last if $chr ne $chrc;

        my $coor = $chr.':'.$pos;
        if ($coor eq $prevCoor) {
          $djindex += 1;
        } else {
          $djindex = 1;
        }

        if ( $vard > 0 and $depth > 0 ) {
          if ( $OR{$coor}{$djindex} ne '' ) {
            my @information = split(",", $OR{$coor}{$djindex});
            my $alt = $information[2];
            my $altd;
            if ($alt eq 'A') {
              $altd = $A;
            } elsif ($alt eq 'C') {
              $altd = $C;
            } elsif ($alt eq 'G') {
              $altd = $G;
            } elsif ($alt eq 'T') {
              $altd = $T;
            } else {
              print STDERR "$coor\t$alt\talt is not ACGT\n";
              exit 22;
            }
            if ($altd > 0) {
              if (exists($blood{$name}) or $blood eq 'yes'){ #it is blood
                $somatic{$coor}{$djindex}{$name} = sprintf("%.3f", $altd/$depth);
              } else {  #it is tumor
                my $endratio = sprintf("%.4f", $vends/$vard);
                if (($endratio <= 0.8 or ($altd - $vends) >= 2) and (($cmean < 3 and $cmedian <= 3) or ($cmean <= 3 and $cmedian < 3))) {  #limiting endsratio and mismatch stuff
                  $somatic{$coor}{$djindex}{$name} = sprintf("%.3f", $altd/$depth);
                } else {  #looks like artifact
                  $somatic{$coor}{$djindex}{$name} = 0;
                  $cmean = 0; #reset for artifact like stuff
                  $cmedian = 0; #reset
                }
              }
            } else {
              $somatic{$coor}{$djindex}{$name} = 0;
            }
          } else {   #coor not found
            if (exists($blood{$name}) or $blood eq 'yes'){ #it is blood
              $somatic{$coor}{$djindex}{$name} = sprintf("%.3f", max($A,$C,$G,$T)/$depth);
            } else { #it is tumor
              my $endratio = sprintf("%.4f", $vends/$vard);
              if (($endratio <= 0.8 or ($vard - $vends) >= 2) and (($cmean < 3 and $cmedian <= 3) or ($cmean <= 3 and $cmedian < 3))) { #limiting endsratio and mismatch stuff
                $somatic{$coor}{$djindex}{$name} = sprintf("%.3f", max($A,$C,$G,$T)/$depth);
              } else {          #looks like artifact
                $somatic{$coor}{$djindex}{$name} = 0;
                $cmean = 0;     #reset for artifact like stuff
                $cmedian = 0;   #reset
              }
            }
          }
        } else {
          $somatic{$coor}{$djindex}{$name} = 0;
        }
        $somatic{$coor}{$djindex}{$name} .= "\t$depth";
        if ($junction != 0) {    #there are some junction reads
          $somatic{$coor}{$djindex}{$name} .= ",$junction";
        }
        $somatic{$coor}{$djindex}{'info'} = join("\t", ("ref","alt"));
        $somatic{$coor}{$djindex}{'consecutive'} .= $cmean.','.$cmedian.','.$vard.';';

        #redefine
        $prevCoor = $coor;
      } #snv
    }
    close IN;
  }

  foreach my $coor (sort {$a =~ /^(\w+):(\d+)$/; my $ca = $1; my $pa = $2; $b =~ /^(\w+):(\d+)$/; my $cb = $1; my $pb = $2; $ca cmp $cb || $pa <=> $pb} keys %somatic) {
    $coor =~ /^(\w+):(\d+)$/;
    my $chrom = $1;
    my $pos = $2;

    foreach my $djindex (sort {$a <=> $b} keys %{$somatic{$coor}}) { #each identical coordinates
      my $info = $somatic{$coor}{$djindex}{'info'};
      my @consecutive = split (';', $somatic{$coor}{$djindex}{'consecutive'});
      my $n = 0;
      my $sumCmean = 0;
      my $sumCmedian = 0;

      foreach my $consecutive (@consecutive) {
        next if $consecutive eq '';
        my @tmp = split (',', $consecutive);
        next if $tmp[0] == 0;
        next if $tmp[1] == 0;
        if ( $tmp[2] >= 3 ) { #vard is the third value of the array, at least three vard to give weights
          $sumCmean += $tmp[0];
          $sumCmedian += $tmp[1];
          $n++;
        }
      }

      if ($n > 0) {         #if you have cmean and cmedian information
        $sumCmean = sprintf("%.2f", ($sumCmean/$n));
        $sumCmedian = sprintf("%.2f", ($sumCmedian/$n));
      }

      my @information = split(",", $OR{$coor}{$djindex});
      my $id = $information[0];
      if ($type =~ 'snv') {
        $info = $information[1]."\t".$information[2];
      }
      next if ($id eq '');

      print "$chrom\t$pos\t$id\t$info";
      foreach my $name (sort {$a =~ /$prefix(\d+)/; my $ia = $1; $b =~ /$prefix(\d+)/; my $ib = $1; $ia <=> $ib} keys %samples) {
        if ($somatic{$coor}{$djindex}{$name} ne '') {
          print "\t$somatic{$coor}{$djindex}{$name}";
        } elsif ($blood eq 'yes') {
          print "\t0\t0";
        }
      }
      print "\t$sumCmean\t$sumCmedian";
      print "\n";
    } #each identical coordinates
  } #each coor
} #each chr

exit 0;



sub getchrpos {

  my $dbfile = shift;

  my $chr_old = "UNDEF";

  my %chr_start;

  my $jumper = 0;

  open DBFILE, "$dbfile" or die "The db file read error!";

  while ( <DBFILE> ) {

    if ($_ =~ /^[\#\@]/){
        $jumper = tell DBFILE;
        next;
    }
    chomp;

    my @cols = split /\t/;
    my $chr  = $cols[0];

    #$chr = 'chr'.$chr unless ($chr =~ /^chr/);
    #$chr = 'chrM' if ($chr eq 'chrMT');

    if ($chr ne $chr_old) {
        $chr_start{$chr} = $jumper;
    }

    $chr_old = $chr;
    $jumper = tell DBFILE;
  }

  close DBFILE;
  print STDERR "$dbfile chr_start loaded\n";

  return \%chr_start;
}
