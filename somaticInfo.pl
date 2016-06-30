use strict;
use Data::Dumper;
use File::Glob ':glob';
use File::Basename;

my $dir = shift;
my @files = bsd_glob("$dir/*.bam");
my %data;

foreach my $bam (@files) {
  $bam = basename($bam);
  $bam =~ /([A-Z0-9\-]+)\_([A-Za-z0-9]+)\.bam/;
  my $id = $1;
  my $type = $2;
  if ($type =~ /Tumor/) {
    push(@{$data{$id}{'T'}}, $bam)
  } else {
    push(@{$data{$id}{'N'}}, $bam)
  }
}

#print STDERR Dumper(\%data);

foreach my $id (keys %data) {
  foreach my $tumor (@{$data{$id}{'T'}}) {
    if ($data{$id}{'N'} ne '') {
      my @normal = @{$data{$id}{'N'}};
      if ($#normal > 0) {
        die("more than one normal for $id");
      }
      print "$tumor\t$normal[0]\n";
    } else {
      print STDERR "$id does not have normal!!!\n";
    }
  }
}

exit 0;
