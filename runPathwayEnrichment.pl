#!/usr/bin/perl -w

use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;
use FindBin qw($RealBin);

my $file = shift;
my $pathway = shift;
my $path = dirname($file);
my $bin = "/tools/trick/";

my $cmd = "perl $bin/pathwayEnrichment.pl $file $pathway >$file\.path";
RunCommand($cmd, 0, 0);

$cmd = "cut -f 6,7,8 $file\.path >ptmp";
RunCommand($cmd, 0, 0);

$cmd = "R CMD BATCH --no-save --no-restore "."\'--args path=\"$path\" file=\"$path\/ptmp\"' $bin/phyper.R $path/R\_html\.out";
RunCommand($cmd, 0, 0);

$cmd = "paste $file\.path $path/pvalues | sort -t \"\t\" -k 10,10 -g >$file\.path\.result";
RunCommand($cmd, 0, 0);

$cmd = "rm $path/ptmp $path/pvalues -f";
RunCommand($cmd, 0, 0);

sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}
