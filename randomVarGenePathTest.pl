use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;
use FindBin qw($RealBin);

my $bin = $RealBin;

my $vcf = shift;       #this is the dbSNP rare variant database
my $pathres = shift;
my $outdir = shift;
my $times = shift;
my $pathdb = shift;
my $Rbin = shift;
my $nonexonic = shift;

if ($vcf eq '' or $vcf =~ /\-h(elp)?$/){
  print STDERR "arguments: rareVardbSNP  pathres  outdir  times  pathdb  Rbin <nonexonic>\n";
  exit 22;
}

open IN, "$vcf";
my $i = 1;
my %vcf;
while ( <IN> ) {
  chomp;
  if (/^#/) {
    next;
  } else {
    if ($nonexonic eq '') {
      $vcf{$i} = $_;
    } else {
      my @cols = split /\t/;
      $vcf{$i} = $cols[1];
    }
    $i++;
  }
}
close IN;
print STDERR "$vcf loaded\n";


my %result;
my $querysize = 0;
if (-e "$pathres") {
  open P, "$pathres";
  while ( <P> ) {
    chomp;
    my @cols = split /\t/;
    my $pathway = join("\t", $cols[0],$cols[1],$cols[2],$cols[3]);
    $result{$pathway}{'genes'} = $cols[4];
    $result{$pathway}{'size'} = $cols[5];
    $result{$pathway}{'mapsize'} = $cols[7];
    $result{$pathway}{'qsize'} = $cols[6];
    $querysize = $result{$pathway}{'qsize'} if ($querysize == 0);
    $result{$pathway}{'better'} = 0;
  }
  close P;
  print STDERR "$pathres loaded\n";
}


srand();

my $totalG = $querysize;
my $totalV = $i - 1;
my %randon;
for (1..$times) {  #randomization 5000 times

  print STDERR "$_\t$_\t$_\t$_\n";
  my $gcount;
  my %lines;
  my %geneCount;
  open OUT, ">$outdir/rgenelist";
  while (1) {  #produce gene list
    my $line = int(rand($totalV));
    while ( exists($lines{$line}) or $line == 0 ) {
      $line = int(rand($totalV));
    }
    $lines{$line} = '';  #remember this line
    #my @genes = &grepGene($vcf{$line});
    #$gcount += scalar(@genes);
    my $genenow;
    if ($nonexonic eq '') {
      my @cols = split(/\t/, $vcf{$line});
      $cols[2] =~ /^(.+?)\:/;
      $genenow = $1;
    } else {
      my @genesnow = split(/\,/, $vcf{$line});
      $genenow = $genesnow[0];
    }
    if (!exists($geneCount{$genenow})){
      $geneCount{$genenow} = '';
      $gcount ++;
    } else {
      next;
    }
    if ($gcount <= $totalG) {
      #foreach my $genenow (@genes) {
      print OUT "$genenow\n";
      #}
    } else {
      last;
    }
  }
  close OUT.
  print STDERR "$outdir/rgenelist is generated.\n";

  #do enrichment analysis
  my $cmd = "perl $bin/runPathwayEnrichment.pl $outdir/rgenelist $pathdb $Rbin";
  RunCommand($cmd, 0 , 0);

  #save it into memory
  open PATH, "$outdir/rgenelist.path.result";
  while ( <PATH> ) {
    chomp;
    my @cols = split /\t/;
    my $pathway = join("\t", $cols[0],$cols[1],$cols[2],$cols[3]);
    if ($cols[7] > $result{$pathway}{mapsize}) {
      $result{$pathway}{'better'} ++;
    }
  }
  close PATH;

  $cmd = "rm $outdir/rgenelist.path $outdir/rgenelist.path.result -f";
  RunCommand($cmd, 0 , 0);
}

foreach my $pathway (sort {$result{$a}{'better'} <=> $result{$b}{'better'} } keys %result){
  my $p = sprintf("%.5f", $result{$pathway}{'better'}/5000);
  print "$pathway";
  print "\t$result{$pathway}{'genes'}";
  print "\t$result{$pathway}{'size'}";
  print "\t$result{$pathway}{'qsize'}";
  print "\t$result{$pathway}{'mapsize'}";
  print "\t$result{$pathway}{'better'}";
  print "\t$p\n";
}

exit 0;


sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}


sub grepGene {
  my $line = shift;
  my @dummy;
  if ($line =~ /function=exonic\;functionalClass=([^\;]+)\;geneName=([\w\.\-\,\;\/]+?)\;[\w\.]+\=/) {
    my $g1 = $2;
    my $function = $1;
    if ($function =~ /^synonymous/ or $function =~ /^nonframeshift/) {
      return(@dummy);
    } else {
      &splitGene($g1);
    }
  } else {
    print STDERR "error:\t\t$_\n";
  }
}


sub splitGene {
  my $genes = shift;
  $genes =~ s/\(.+?\)//g;                #remove the splicing parentheses
  my @genes = split (/[\;\,]/, $genes);
  my %tmp;
  foreach my $gene (@genes) {
    $tmp{$gene} = "";
  }
  @genes = keys %tmp;
  return(@genes);
}
