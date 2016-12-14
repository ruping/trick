use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;
use FindBin qw($RealBin);

my $bin = $RealBin;

my $files = shift;
my @files = split(',', $files);

my $pathdb = shift;
my $Rbin = shift;



foreach my $file (@files) {

  my $cmd = "perl $bin/runPathwayEnrichment.pl $file $pathdb $Rbin";
  RunCommand($cmd, 0 , 0);

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

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}
