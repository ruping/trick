use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;
use FindBin qw($RealBin);

my $bin = $RealBin;

my $rare = shift;
my $samples = shift;
my $outdir = shift;
my $pathdb = shift;
my $Rbin = shift;
my @samples = split(',', $samples);

foreach my $sample (@samples) {
  my $lastcolindex = `perl $bin/columnIndex.pl $sample $rare`;
  my $cindex =~ s/\n$//;
  $cindex += 1;
  my $cmd = "awk \-F\"\\t\" \'$cindex != 0\' $rare >$outdir/$sample\.rare";
  RunCommand($cmd, 0, 0);

  #do enrichment analysis
  $cmd = "perl $bin/runPathwayEnrichment.pl $outdir/$sample\.rare $pathdb $Rbin";
  RunCommand($cmd, 0, 0);
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
