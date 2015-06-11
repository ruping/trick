use strict;
use File::Glob ':glob';
use Data::Dumper;
use Getopt::Long;
use FindBin qw($RealBin);

my $list;   #filename of all vcfs
my $type;
my $normal;
my $recheck;
my $task;
my $prefix;
my $dbsnp = "no";
my $clinical;
my $tmpdir = "./";
my $bin = $RealBin;
my $tolerance = 0;
my $strandBiasTh = 0.005;
my $tailDisBiasTh = 0.005;

GetOptions (
            "list|l=s"       => \$list,             #filename of all vcfs
            "type|t=s"       => \$type,             #snv or indel
            "normal|n=s"     => \$normal,           #comma seperated id of normal samples
            "recheck|r=s"    => \$recheck,          #directory of rechecked files
            "task|k=s"       => \$task,             #task type
            "prefix|p=s"     => \$prefix,
            "dbsnp|d=s"      => \$dbsnp,
            "clinical|c=s"   => \$clinical,         #clinical dbSNP sites
            "tmpdir|y=s"     => \$tmpdir,
            "tolerance=i"    => \$tolerance,
            "strandBiasTh=f" => \$strandBiasTh,
            "tailDisBiasTh=f"=> \$tailDisBiasTh,
            "help|h"         => sub{
                               print "usage: $0 get all somatic and rare variants from a bunch of vcf files\n\nOptions:\n\t--list\t\tthe filename of all vcfs\n";
                               print "\t--type\t\tthe type of variants, snv or indel\n";
                               print "\t--normal\tcomma seperated id of normal samples\n";
                               print "\t--prefix\tthe prefix of samples' names\n";
                               print "\t--task\t\tthe task, such as tcga or rnaediting\n";
                               print "\t--recheck\tthe dir whether recheck files are located\n";
                               print "\t--dbsnp\t\tyes or no, whether to keep dbsnp variants into the table\n";
                               print "\t--clinical\tthe file containing the clinical dbSNP sites, it is a gzipped vcf file\n";
                               print "\t--tolerance\tthe tolerance for comparing indels with certain distance shift\n";
                               print "\t--tmpdir\tthe temporary dir to write tmp files\n";
                               print "\t--strandBiasTh\tthe p value threshold for filtering out vars with strand bias <0.005>\n";
                               print "\t--tailDisBiasTh\tthe p value threshold for filtering out vars with Tail distance bias <0.005>\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );


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
  if ($task ne 'tcga' and $file =~ /($prefix\d+)[^a-zA-Z0-9]/) {
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

  if ($clinical ne ''){        # want to grep the clinical sites only, first generate tmp intersect file
    my $cmd = "perl $bin/intersectFiles.pl -o $file -m $clinical -t $tolerance -vcf -overlap -column INFO >$tmpdir/tmp";
    RunCommand($cmd,0,0);
    $file = "$tmpdir/tmp";
  }


  my $openway = $file;
  if ($file =~ /\.gz$/){
    $openway = "gzip -dc $file |";
  }
  print STDERR "$openway\n";

  open IN, "$openway";
  my $revertornot = "no";
  my $printerror = 0;
  my $singlecalling = "no";
  while ( <IN> ) {
     chomp;
     if ($_ =~ /^#/) {
       if ($_ =~ /^#CHROM\tPOS\tID/) {
         my @cols = split /\t/;
         my $minusI = ($clinical eq '')? 1:2;
         if ($cols[$#cols - $minusI] eq $normal) {
           $revertornot = "yes";
         } elsif ( $cols[$#cols - $minusI] eq 'FORMAT' ) {
           $singlecalling = "yes";
         }
         print STDERR "revert or not: $revertornot\n";
         print STDERR "singlecalling: $singlecalling\n";
         next;
       } else {
         next;
       }
     }
     my ($chr, $pos, $id, $ref, $alt, $qual, $pass, $info, $format, $sample, $blood, $clinINFO) = split /\t/;

     next if ($qual ne '.' and $qual < 30);

     if ($revertornot eq 'yes') {   #revert sample and blood
        my $tmp = $sample;
        $sample = $blood;
        $blood = $tmp;
     }

     my @formats = split(':', $format);
     my %formindex;
     for(my $f = 0; $f <= $#formats; $f++) {
       $formindex{$formats[$f]} = $f;
     }
     if ($printerror == 0){
       print STDERR Dumper(\%formindex);
       $printerror ++;
     }

     my @sample = split(/\:/,$sample);

     if ($sample[$formindex{'GT'}] !~ /1/) {     #skip some wierd thing
        next;
     }

     if ($clinical ne ''){  # want to grep the clinical sites only
        goto PRODUCE;
     }

     ###########################################################################decide somatic
     my $somatic = 0;
     if ($type eq 'snv') {   #for snp
        my @blood = split(/\:/,$blood);
        my @bad = split (/\,/, $blood[$formindex{'AD'}]);
        if (scalar(@bad) > 1) {
           if ($blood[$formindex{'GT'}] !~ /1/ and $bad[1] == 0) {
             $somatic = 1;
           }
        } else {
           if ($blood[$formindex{'GT'}] !~ /1/) {
             $somatic = 1;
           }
        }
        if ($recheck ne '') {       #recheck
           my $coor = $chr.':'.$pos;
           if (exists($recheck{$coor})) {
              $somatic = 1;
           }
        }
     } elsif ($type eq 'indel') {
         my @blood = split(/\:/,$blood);
         my @bad = split (/\,/, $blood[$formindex{'AD'}]);
         if (scalar(@bad) > 1){
           if ($blood[$formindex{'GT'}] !~ /1/ and $bad[1] == 0) {
             $somatic = 1;
           }
         } else {
           if ($blood[$formindex{'GT'}] !~ /1/) {
             $somatic = 1;
           }
         }
         if ($recheck ne '') {        #recheck
           my $coor = $chr.':'.$pos;
           if (exists($recheck{$coor})) {
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
            next if $dbsnp eq "no";
         }
         elsif ($freq > 0 and ($id !~ /^1KG/ and $id !~ /^ESP5400/)) {   #still dbSNP ones
            next if $dbsnp eq "no";
         }
       }
     } #somatic common snp

   PRODUCE:

     unless ($qual ne '.' and $qual > 90) { #unless very high quality

       if ($info =~ /MQ0F=(.+?);/) {
         next if ($1 > 0.1);
       }

       if ($info =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/) {
         my $strandb = $1;
         my $baseqb = $2;
         my $mapqb = $3;
         my $tailb = $4;
         if ($strandb =~ /e/) { #strand bias
           next;
         } elsif ($strandb < $strandBiasTh) { #strand bias make it very stringent! for net data
           next;
         } elsif ($baseqb =~ /e/) { #basequality bias
           next;
         } elsif ($baseqb < 0.0005) { #basequality bias
           next;
           #} elsif ($mapqb =~ /e/) {          #mapquality bias   #mask it for sid's data
           #  next;
           #} elsif ($mapqb < 0.0001) {        #mapquality bias   #mask it for sid's data
           #  next;
         } elsif ($tailb =~ /e/) { #tailbias
           next;
         } elsif ($tailb < $tailDisBiasTh) { #tailbias
           next;
         } else {
           #pass
         }
       }

     }  #unless very high quality

     my $coor = $chr.':'.$pos;

     my $maf = -1;
     if ($info =~ /DP4\=(\d+)\,(\d+)\,(\d+)\,(\d+)\;/) {
        $maf = sprintf("%.3f", ($3+$4)/($1+$2+$3+$4));
     } else {
        my @sampleinfo = split(":", $sample);
        $sampleinfo[$formindex{'AD'}] =~ /^(\d+)\,(\d+)$/;
        $maf = sprintf("%.3f", $2/($1+$2));
     }

     $somatic{$coor}{$name} = $maf.'|'.$qual;
     my $function;
     $info =~ /(function=.+?$)/;
     $function = $1;
     #$somatic{$coor}{'function'} = $function;
     my $idrefalt = join("\t", ($id,$ref,$alt));
     $somatic{$coor}{'info'}{$idrefalt} = $function;  #not necessarily just one!

     if ($somatic == 1) {
        $somatic{$coor}{'somatic'} .= $name.',';
     } else {
        $somatic{$coor}{'germline'} .= $name.',';
     }
     if ($clinical ne '') {
       $somatic{$coor}{'clinical'} = $clinINFO;  #clinical id information
     }
  }
  close IN;
  $tmpc = scalar(keys %somatic);
  print STDERR "table tmp: $tmpc\n";

  my $cmd = "rm $tmpdir/tmp -f";
  RunCommand($cmd, 0 ,0);
}

#######------------------print header-----------------------#####
print "#chr\tpos\tid\tref\talt";
foreach my $name (sort {$a =~ /$prefix(\d+)/; my $ia = $1; $b =~ /$prefix(\d+)/; my $ib = $1; $ia <=> $ib} keys %samples) {
  print "\t$name";
}
print "\tfunction";
if ($clinical eq '') {
  print "\tsomatic\tgermline";
} else {
  print "\tclinical";
}
print "\n";
#################################################################


foreach my $coor (sort {$a =~ /^(\w+):(\d+)$/; my $ca = $1; my $pa = $2; $b =~ /^(\w+):(\d+)$/; my $cb = $1; my $pb = $2; $ca cmp $cb || $pa <=> $pb} keys %somatic) {
  $coor =~ /^(\w+):(\d+)$/;
  my $chrom = $1;
  my $position = $2;
  foreach my $info (keys (%{$somatic{$coor}{'info'}})) {
    print "$chrom\t$position\t$info";
    foreach my $name (sort {$a =~ /$prefix(\d+)/; my $ia = $1; $b =~ /$prefix(\d+)/; my $ib = $1; $ia <=> $ib} keys %samples) {
      if ($somatic{$coor}{$name} ne '') {
        print "\t$somatic{$coor}{$name}";
      } else {
        print "\t0";
      }
    }
    my $function = $somatic{$coor}{'info'}{$info};
    my $somatic = ($somatic{$coor}{'somatic'} eq '')? 0 : $somatic{$coor}{'somatic'};
    my $germline = ($somatic{$coor}{'germline'} eq '')? 0 : $somatic{$coor}{'germline'};
    print "\t$function";
    if ($clinical eq '') {
      print "\t$somatic\t$germline";
    } else {
      print "\t$somatic{$coor}{'clinical'}";
    }
    print "\n";
  } #for each different variant in this coordinate
}


sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}
