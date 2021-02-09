use strict;
use File::Glob ':glob';
use Data::Dumper;
use File::Basename;

my $dir = shift;
my @tfths = bsd_glob("$dir/*_thresholds.txt");
#print STDERR Dumper(\@tfths);

foreach my $tfth (@tfths){
  $tfth =~ /$dir[\/]+(.+?)\_thresholds.txt/;
  my $tf = $1;
  #print STDERR "$tf\n";
  my $cmd = "R CMD BATCH --no-save --no-restore "."\'--args path=\"$dir\" tf=\"$tf\"\' /tools/trick/TFthresLinear.R R_tmp.out";
  system($cmd);
}

exit 0;
