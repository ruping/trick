use strict;
use File::Glob ':glob';


my $dir = shift;
my $suffix = shift;

my @files = bsd_glob("$dir/*");

foreach my $file (@files){
  my $newfile = $file;
  $newfile .= "\.$suffix";
  my $cmd = "mv $file $newfile";
  system($cmd);
}
