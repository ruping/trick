use strict;
use Data::Dumper;
use File::Glob ':glob';
use File::Basename;

my $dir = shift;
my $suffix = shift;
my $sdir = shift;
my $sdirsuffix = shift;

if ($dir eq '' or $suffix eq '' or $sdir eq '') {
  print "usage: moveRunLogs.pl dir suffix outdir outdirsuffix\n\n";
  exit 22;
}

my @logs = bsd_glob("$dir/*$suffix");
print STDERR Dumper(\@logs);
print STDERR "$sdir\t$sdirsuffix\n";

foreach my $log (@logs) {
  my $b = basename($log);
  $b =~ /^(.+?)\.$suffix$/;
  my $sample = $1;

  print STDERR "$sample\n";

  my $ddir = $sdir;
  $ddir .= $sample.'/'.$sdirsuffix;
  my $sfile = $ddir.$b;

  print STDERR "original file: $log\ndestination file: $sfile\n";

  my $cmd = "mv $log $sfile";

  unless (-e "$sfile") {
    system($cmd);
  }
}

exit 0;
