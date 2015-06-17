use strict;
use File::Glob ':glob';
use File::Basename;

my $dir = shift;
my $pattern = shift;
my $changed = shift;
my @files = bsd_glob("$dir/*$pattern*");

foreach my $file (@files){
  my $basename = basename($file);
  my $dirname = dirname($file);
  #(my $newfile = $file) =~ s/$pattern/$changed/;
  $basename =~ s/$pattern/$changed/;
  my $newfile = $dirname.'/'.$basename;
  my $cmd = "mv -f $file $newfile";
  system($cmd);
}
