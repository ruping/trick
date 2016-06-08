use strict;
use File::Glob ':glob';
use File::Basename;

my $dir = shift;
my $pattern = shift;
my $outdir = shift;

my @files = bsd_glob("$dir/*");

foreach my $f (@files){
  my $b = basename($b);
  next if $b !~ /$pattern/;
  my $cmd = "touch $outdir/$b";
  system($cmd);
}

exit 0;
