use strict;
use Data::Dumper;
use File::Glob ':glob';
use File::Basename;

my $dir = shift;
my $suffix = shift;
my $sdir = shift;
my $sdirsuffix = shift;

if ($dir eq '' or $suffix eq '' or $sdir eq ''){
  print "usage: moveRunLogs.pl dir suffix outdir outdirsuffix\n\n";
  exit 22;
}

my @logs = bsd_glob("$dir/*$suffix");

foreach my $log (@logs) {
  my $b = basename($log);
  $b =~ /^(.+?)\.$suffix$/;
  my $sample = $1;

  $sdir .= $sample.'/'.$sdirsuffix;
  my $sfile = $sdir.$b;
  my $cmd = "mv $log $sfile";
  if (-e "$sfile"){
    system($cmd);
  }
}

exit 0;
