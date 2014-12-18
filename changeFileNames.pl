use strict;
use File::Glob ':glob';

my $dir = shift;
my $pattern = shift;
my $changed = shift;
my @files = bsd_glob("$dir/*$pattern*");


foreach my $file (@files){
  (my $newfile = $file) =~ s/$pattern/$changed/;
  my $cmd = "mv -f $file $newfile";
  system($cmd);
}
