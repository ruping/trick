use strict;
use File::Glob ':glob';
use File::Basename;
use Data::Dumper; 

my $dir = shift;
my $pattern = shift;
my $changed = shift;

my @files = bsd_glob("$dir/*");
print STDERR "pattern: $pattern\n";
print STDERR "changed: $changed\n";


foreach my $file (@files){
  next if $files !~ /$pattern/;
  my $basename = basename($file);
  my $dirname = dirname($file);
  $basename =~ s/$pattern/$changed/;
  my $newfile = $dirname.'/'.$basename;
  my $cmd = "mv -f $file $newfile";
  system($cmd);
}
