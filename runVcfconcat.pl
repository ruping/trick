use strict;
use File::Glob ':glob';
use File::Basename;

my $dir = shift;
my @vcfs = bsd_glob("$dir/*.vcf");

my %vcfs;
foreach my $vcf (@vcfs){
  my $basename = basename($vcf);
  if ($basename =~ /^([A-Za-z0-9\-\_]+)\.([A-Za-z0-9]+)\./){
    my $sample = $1;
    my $chr = $2;
    $vcfs{$chr} = $vcf;
  } else {
    print "strang $vcf\n";
    exit 22;
  }
}

my $vcfs = '';
foreach my $chr (sort {$a cmp $b} keys %vcfs){
  $vcfs .= "$vcfs{$chr} ";
}

print "$vcfs\n";
