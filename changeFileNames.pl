use strict;
use File::Glob ':glob';
use File::Basename;
use Data::Dumper; 

my $dir = shift;
my $pattern = shift;
my $changed = shift;

my @files = bsd_glob("$dir/*$pattern*");
print STDERR "pattern: $pattern\n";
print STDERR "changed: $changed\n";
print STDERR (\@files);

foreach my $file (@files){
  my $basename = basename($file);
  my $dirname = dirname($file);
  $basename =~ s/$pattern/$changed/;
  my $newfile = $dirname.'/'.$basename;
  my $cmd = "mv -f $file $newfile";
  system($cmd);
}
