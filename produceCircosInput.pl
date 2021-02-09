use strict;
use File::Glob ':glob';
use File::Basename;

my $dir = shift;
my $outdir = shift;

my @segments = bsd_glob("$dir/*_nclones1.TitanCNA.segments.txt");
foreach my $seg (@segments){
  my $b = basename($seg);
  $b =~ /^([A-Z0-9a-z\_\-]+)_nclones1\.TitanCNA\.segments\.txt/;
  my $samp = $1;
  if ($samp eq ''){
    print STDERR "$b is wierd naming!\n";
    exit 22;
  }

  my $cmd = "perl /tools/trick/prepareTitanCircos.pl $seg >$outdir/$samp\_nt 2>$outdir/$samp\_nb";
  system($cmd);
  print STDERR "$b finished\n";
}

exit 0;
