use strict;
use File::Glob ':glob';
use Data::Dumper;

my $dir = shift;
my $pattern = shift;
my $dens = shift;
if ($dens eq ''){
  $dens = 300;
}

my @files = bsd_glob("$dir/*$pattern*.pdf");

foreach my $file (@files) {

  my $out = $file;
  $out =~ s/\.pdf$/\.png/;

  my $cmd = "convert -density $dens $file -background White -flatten -quality 95 $out";
  print STDERR "$cmd\n";
  system($cmd);
}

exit 0;
