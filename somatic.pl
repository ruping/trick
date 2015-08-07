use strict;
use File::Glob ':glob';
use Data::Dumper;
use Getopt::Long;
use FindBin qw($RealBin);
use File::Basename;
use List::Util qw(max min sum);

my $list;   #filename of all vcfs
my $type;
my $normal;
my $task;
my $prefix;
my $dbsnp = "no";
my $tmpdir = "./";
my $bin = $RealBin;
my $tolerance = 0;
my $strandBiasTh = 0.005;
my $tailDisBiasTh = 0.005;
my $nonsegdup;
my $exonic;

GetOptions (
            "list|l=s"       => \$list,             #filename of all vcfs
            "type|t=s"       => \$type,             #snv or indel
            "normal|n=s"     => \$normal,           #comma seperated id of normal samples
            "task|k=s"       => \$task,             #task type
            "prefix|p=s"     => \$prefix,
            "dbsnp|d=s"      => \$dbsnp,
            "tmpdir|y=s"     => \$tmpdir,
            "tolerance=i"    => \$tolerance,
            "strandBiasTh=f" => \$strandBiasTh,
            "tailDisBiasTh=f"=> \$tailDisBiasTh,
            "nonsegdup"      => \$nonsegdup,
            "exonic"         => \$exonic,
            "help|h"         => sub{
                               print "usage: $0 get all somatic and rare variants from a bunch of vcf files\n\nOptions:\n\t--list\t\tthe filename of all vcfs\n";
                               print "\t--type\t\tthe type of variants, snv or indel\n";
                               print "\t--normal\tcomma seperated id of normal samples\n";
                               print "\t--prefix\tthe prefix of samples' names\n";
                               print "\t--task\t\tthe task, such as tcga or rnaediting\n";
                               print "\t--dbsnp\t\tyes or no, whether to keep dbsnp variants into the table\n";
                               print "\t--tolerance\tthe tolerance for comparing indels with certain distance shift\n";
                               print "\t--tmpdir\tthe temporary dir to write tmp files\n";
                               print "\t--strandBiasTh\tthe p value threshold for filtering out vars with strand bias <0.005>\n";
                               print "\t--tailDisBiasTh\tthe p value threshold for filtering out vars with Tail distance bias <0.005>\n";
                               print "\t--nonsegdup\tdo not collect segdup ones\n";
                               print "\t--exonic\tonly collect exonic ones (including UTR and splicing)\n";
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

my @prefix = split(',', $prefix);
my $prefixReg = join('|', @prefix);
print STDERR "prefixReg is $prefixReg\n";

my %somatic;
my %samples;
foreach my $file (@list) {
  my $name;
  my $filebase = basename($file);
  if ($prefixReg ne '' and $filebase =~ /(($prefixReg)[a-zA-Z0-9\-\_]+)[^a-zA-Z0-9\-\_]/) {     #changed to adapt to lowGI data
    $name = $1;
  }
  elsif ($task =~ /tcga/i) {
    #$file =~ /\/((TCGA\-[^\-]+\-[^\-]+)\-[^\-]+\-[^\-]+\-[^\-]+\-\d+)\./;
    $filebase =~ /(TCGA\-[A-Z0-9]+\-[A-Z0-9]+)/;
    $name = $1;
  }
  print STDERR "$name\t$file\t$filebase\n";

  $samples{$name} = '';

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
       if ($_ =~ /^#CHROM\tPOS\tID/) {   #the common three column header in vcf file
         my @cols = split /\t/;
         my $minusI = 1;
         if (($cols[$#cols - $minusI] eq $normal) or ($cols[$#cols - $minusI] =~ /NORMAL/i)) {
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
     my ($chr, $pos, $id, $ref, $alt, $qual, $pass, $info, $format, $sample, $blood) = split /\t/;

     next if ($qual ne '.' and $qual < 30 and $pass ne 'PASS');

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
     if ($printerror == 0) {
       print STDERR Dumper(\%formindex);
       $printerror ++;
     }

     my @sample = split(/\:/,$sample);

     if ($sample[$formindex{'GT'}] !~ /1/) {     #skip some wierd thing
        next;
     }


     ###########################################################################decide somatic
     my $somatic = 0;
     if ($singlecalling eq 'no'){      # if it is paired calling
       if ($type eq 'snv') {    #for snp
         my @blood = split(/\:/,$blood);
         if ($blood[$formindex{'GT'}] !~ /1/) {
           $somatic = 1;
         }
       } elsif ($type eq 'indel') { #for indel
         my @blood = split(/\:/,$blood);
         if ($blood[$formindex{'GT'}] !~ /1/) {
           $somatic = 1;
         }
       }
     }
     #############################################################################decide somatic

     #print STDERR "somatic: $somatic\n";

     if ($id ne '.' or $info =~ /dbSNP/ or $info =~ /1KG\=/ or $info =~ /ESP\d+\=/) {  #snp in population, is it a somatic one?

       my $freq = -1;
       if ($info =~ /(1KG=(.+?));/) {
           my $kid = $1;                               #re define $id to 1KG when absent
           $freq = $2;
           if ($id eq '.') {
              $id = $kid;
           }
       }

       if ($info =~ /((ESP\d+)=(.+?));/) {
           my $kid = $1;                               #re define $id to ESP when absent
           $freq = $3;
           if ($id eq '.') {
              $id = $kid;
           }
       }

       if ( $somatic == 0 ) {                            #if it is not somatic, then only rare ones should be kept
         if ($freq == -1) {
           if ($dbsnp eq "no" or $task =~ /rare/i) {
             next;
           }
         } else {   #freq is defined
           if ($id !~ /^1KG/ and $id !~ /^ESP\d+/) {  #dbSNP ones with reported MAF in 1KG or ESP
             next if $dbsnp eq "no";
           }
           if ($task =~ /rare/i) {
             next if $freq > 0.005;                   #rare kept
           }
         }
       }
     } #somatic common snp

   PRODUCE:

     unless ($qual ne '.' and $qual > 90) { #unless very high quality

       if ($info =~ /MQ0F=(.+?);/) {    #if with MQ0F
         next if ($1 > 0.1);
       }

       if ($info =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/) {    #if with PV4 info
         my $strandb = $1;
         my $baseqb = $2;
         my $mapqb = $3;
         my $tailb = $4;
         if ($strandb =~ /e/) {               #strand bias
           next;
         } elsif ($strandb < $strandBiasTh) { #strand bias make it very stringent! for net data
           next;
         } elsif ($baseqb =~ /e/) {           #basequality bias
           next;
         } elsif ($baseqb < 0.0005) {         #basequality bias
           next;
           #} elsif ($mapqb =~ /e/) {         #mapquality bias   #mask it for sid's data
           #  next;
           #} elsif ($mapqb < 0.0001) {       #mapquality bias   #mask it for sid's data
           #  next;/MQ0
         } elsif ($tailb =~ /e/) {            #tailbias
           next;
         } elsif ($tailb < $tailDisBiasTh) {  #tailbias
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
        if ($sampleinfo[$formindex{'AD'}] =~ /^(\d+)\,(\d+)$/){
          $maf = sprintf("%.3f", $2/($1+$2));
        } else {
          my @ads = split(',', $sampleinfo[$formindex{'AD'}]);
          my $adsum = sum(@ads);
          shift(@ads);
          my $altadsum = sum(@ads);
          $maf = sprintf("%.3f", $altadsum/$adsum);
        }
     }

     $somatic{$coor}{$name} = $maf.'|'.$qual;
     my $function;
     $info =~ /(function=.+?$)/;
     $function = $1;

     if ($nonsegdup){
       next if $function =~ /segdup\.score/;
     }
     if ($exonic){
       next if ($function !~ /exonic/ and $function !~ /UTR[35]/ and $function !~ /splicing/);
     }

     my $idrefalt = join("\t", ($id,$ref,$alt));
     $somatic{$coor}{'info'}{$idrefalt} = $function;  #not necessarily just one!

     if ($somatic == 1) {
        $somatic{$coor}{'somatic'} .= $name.',';
     } else {
        $somatic{$coor}{'germline'} .= $name.',';
     }
  }
  close IN;
  my $tmpc = scalar(keys %somatic);
  print STDERR "table tmp: $tmpc\n";

  my $cmd = "rm $tmpdir/tmp -f";
  RunCommand($cmd, 0 ,0);
}

#######------------------print header-----------------------#####
print "#chr\tpos\tid\tref\talt";
if ($prefixReg ne ''){
  foreach my $name (sort {$a =~ /($prefixReg)(\d+)([a-zA-Z0-9\-\_]+)?/; my $ia = $2; my $ias = $3; $b =~ /($prefixReg)(\d+)([a-zA-Z0-9\-\_]+)?/; my $ib = $2; my $ibs = $3; $ia <=> $ib or $ias cmp $ibs} keys %samples) {
    print "\t$name";
  }
} elsif ($task =~ /tcga/i) {
  foreach my $name (sort {$a =~ /TCGA\-([A-Z0-9]+)\-([A-Z0-9]+)/; my $tsa = $1; my $inda = $2; $b =~ /TCGA\-([^\-]+)\-([^\-]+)/; my $tsb = $1; my $indb = $2; $tsa cmp $tsb or $inda cmp $indb} keys %samples){
    print "\t$name";
  }
}
print "\tfunction\ttrace";
print "\n";
#################################################################


foreach my $coor (sort {$a =~ /^(\w+):(\d+)$/; my $ca = $1; my $pa = $2; $b =~ /^(\w+):(\d+)$/; my $cb = $1; my $pb = $2; $ca cmp $cb || $pa <=> $pb} keys %somatic) {
  $coor =~ /^(\w+):(\d+)$/;
  my $chrom = $1;
  my $position = $2;
  foreach my $info (keys (%{$somatic{$coor}{'info'}})) {
    print "$chrom\t$position\t$info";
    if ($prefixReg ne '') {
      foreach my $name (sort {$a =~ /($prefixReg)(\d+)([a-zA-Z0-9\-\_]+)?/; my $ia = $2; my $ias = $3; $b =~ /($prefixReg)(\d+)([a-zA-Z0-9\-\_]+)?/; my $ib = $2; my $ibs = $3; $ia <=> $ib or $ias cmp $ibs} keys %samples) {
        if ($somatic{$coor}{$name} ne '') {
          print "\t$somatic{$coor}{$name}";
        } else {
          print "\t0";
        }
      }
    } elsif ($task =~ /tcga/i) {
      foreach my $name (sort {$a =~ /TCGA\-([A-Z0-9]+)\-([A-Z0-9]+)/; my $tsa = $1; my $inda = $2; $b =~ /TCGA\-([^\-]+)\-([^\-]+)/; my $tsb = $1; my $indb = $2; $tsa cmp $tsb or $inda cmp $indb} keys %samples) {
        if ($somatic{$coor}{$name} ne '') {
          print "\t$somatic{$coor}{$name}";
        } else {
          print "\t0";
        }
      }
    }
    my $function = $somatic{$coor}{'info'}{$info};
    my $somatic = ($somatic{$coor}{'somatic'} eq '')? 0 : $somatic{$coor}{'somatic'};       #temporarily silence somatic and germline info
    my $germline = ($somatic{$coor}{'germline'} eq '')? 0 : $somatic{$coor}{'germline'};    #temporarily silence somatic and germline info
    my $trace = 'somatic='.$somatic.';germline='.$germline;
    print "\t$function\t$trace";
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
