use strict;
use File::Glob ':glob';
use List::Util qw[min max];
use Data::Dumper;

my $file = shift;      #filename of all rechecked files
my $type = shift;
my $original = shift;
my $task = shift;      #only for tcga

my @list;
open IN, "$file";
while ( <IN> ) {
  chomp;
  next if /^#/;
  push(@list, $_);
}
close IN;


####################################################get chr pos
my %chrJumper;
$chrJumper{'original'} = getchrpos($original);

my %samples;
foreach my $file (@list) {
  my $name;
  if ($file =~ /(AC\d+)/) {
    $name = $1;
  } elsif ($task eq 'tcga' and $file =~ /\/((TCGA\-([^\-]+\-[^\-]+))\-[^\-]+\-[^\-]+\-[^\-]+\-\d+)$/) {
    $name = $1;
  }
  $samples{$name} = '';
  $chrJumper{$name} = getchrpos($file);
}

print STDERR Dumper(\%chrJumper);
####################################################get chr pos


print "#chr\tpos\tid\tref\talt";
foreach my $name (sort keys %samples) {
   print "\t$name\t$name".'d';
}
print "\tcmeanav\tcmedianav\n";

foreach my $chrc (sort keys %{$chrJumper{'original'}}) {

  my $jumperO = $chrJumper{'original'}->{$chrc};
  open OR, "$original";
  seek(OR, $jumperO, 0);
  my %OR;
  while ( <OR> ) {
    next if /^[@#]/;
    chomp;
    my @cols = split /\t/;
    my $chr = $cols[0];

    last if $chr ne $chrc;

    my $pos = $cols[1];
    my $id  = $cols[2];
    my $coor = $chr.':'.$pos;
    my $ref = $cols[3];
    my $alt = $cols[4];
    $OR{$coor} = $id.','.$ref.','.$alt;
  }
  close OR;
  print STDERR "$original $chrc loaded\n";


  my %somatic;
  foreach my $file (@list) {
    my $name;
    if ($file =~ /(AC\d+)/) {
      $name = $1;
    } elsif ($task eq 'tcga' and $file =~ /\/((TCGA\-([^\-]+\-[^\-]+))\-[^\-]+\-[^\-]+\-[^\-]+\-\d+)$/) {
      $name = $1;
    }

    print STDERR "$name\n";

    my $jumperI = $chrJumper{$name}->{$chrc};
    open IN, "$file";
    seek(IN, $jumperI, 0);
    while ( <IN> ) {
      chomp;
      next if /^[#@]/;
      next if /^chr\t/;
      if ($type =~ /indel/) {       #indel
        my ($chr, $pos, $ref, $alt, $indelType, $depth, $vard, $junction, $cmean, $cmedian) = split /\t/;
        last if $chr ne $chrc;
        my $coor = $chr.':'.$pos;
        if ($vard > 0 and $depth > 0) {
          $somatic{$coor}{$name} = sprintf("%.3f", $vard/$depth);
        } else {
          $somatic{$coor}{$name} = 0;
        }
        $somatic{$coor}{$name} .= "\t$depth";
        if ($junction != 0) {  #there are some junction reads
          $somatic{$coor}{$name} .= ",$junction";
        }
        $somatic{$coor}{'info'} = join("\t", ($ref,$alt));
        $somatic{$coor}{'consecutive'} .= $cmean.','.$cmedian.';';

      } elsif ($type =~ /snv/) {    #snv

        my ($chr, $pos, $depth, $vard, $A, $C, $G, $T, $vends, $junction, $cmean, $cmedian) = split /\t/;
        last if $chr ne $chrc;
        my $coor = $chr.':'.$pos;
        if ( $vard > 0 and $depth > 0 ) {
          if ( exists( $OR{$coor} ) ) {
            my @information = split(",", $OR{$coor});
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
              $somatic{$coor}{$name} = sprintf("%.3f", $altd/$depth);
            } else {
              $somatic{$coor}{$name} = 0;
            }
          } else {
            $somatic{$coor}{$name} = sprintf("%.3f", max($A,$C,$G,$T)/$depth);
          }
        } else {
          $somatic{$coor}{$name} = 0;
        }
        $somatic{$coor}{$name} .= "\t$depth";
        if ($junction != 0) {    #there are some junction reads
          $somatic{$coor}{$name} .= ",$junction";
        }
        $somatic{$coor}{'info'} = join("\t", ("ref","alt"));
        $somatic{$coor}{'consecutive'} .= $cmean.','.$cmedian.';';
      }  #snv
    }
    close IN;
  }

  foreach my $coor (sort {$a =~ /^(\w+):(\d+)$/; my $ca = $1; my $pa = $2; $b =~ /^(\w+):(\d+)$/; my $cb = $1; my $pb = $2; $ca cmp $cb || $pa <=> $pb} keys %somatic) {
    $coor =~ /^(\w+):(\d+)$/;
    my $chrom = $1;
    my $pos = $2;
    my $info = $somatic{$coor}{'info'};
    my @consecutive = split (';', $somatic{$coor}{'consecutive'});
    my $n = 0;
    my $sumCmean = 0;
    my $sumCmedian = 0;

    foreach my $consecutive (@consecutive) {
       next if $consecutive eq '';
       my @tmp = split (',', $consecutive);
       next if $tmp[0] == 0;
       next if $tmp[1] == 0;
       $sumCmean += $tmp[0];
       $sumCmedian += $tmp[1];
       $n++;
    }

    if ($n > 0) {  #if you have cmean and cmedian information
       $sumCmean = sprintf("%.2f", ($sumCmean/$n));
       $sumCmedian = sprintf("%.2f", ($sumCmedian/$n));
    }

    my @information = split(",", $OR{$coor});
    my $id = $information[0];
    if ($type =~ 'snv') {
       $info = $information[1]."\t".$information[2];
    }
    next if ($id eq '');
    print "$chrom\t$pos\t$id\t$info";
    foreach my $name (sort keys %samples) {
      if ($somatic{$coor}{$name} ne '') {
        print "\t$somatic{$coor}{$name}";
      }
    }
    if ($n > 0) {
       print "\t$sumCmean\t$sumCmedian";
    }
    print "\n";
  }
}

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
