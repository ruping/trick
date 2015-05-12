use strict;
use Data::Dumper;
use File::Glob ':glob';


my $dir = shift;
my $prefix = shift;

my @files = bsd_glob("$dir/$prefix*/FastCallResults_*.txt");
print STDERR Dumper (\@files);


foreach my $file (@files) {

  
}
