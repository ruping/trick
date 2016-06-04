use strict;
use Data::Dumper;
use File::Glob ':glob';
use File::Basename;

my $dir = shift;

my @logs = bsd_glob("$dir/*.run.log");

foreach my $log (@logs){
  my $b = basename($log);
  $b =~ /^(.+?)\.run\.log$/;
  (my $sd = $log) =~ s/\.run\.log$//;
  my $sample = $1;
  my $cmd = "mv $log $sd";
  if (-e "$sd"){
    system($cmd);
  }
}

exit 0;
