use strict;
use File::Glob ':glob';
use Data::Dumper;

my $input = shift;
my $dir = shift;

my @tissue = qw(Rectal_Mucosa Rectal_Smooth_Muscle Duodenum_Mucosa Lung Pancreatic_Islets Colonic_Mucosa);

my $opc = 1;
foreach my $tissue (@tissue) {
  my @files = bsd_glob("$dir/$tissue/*.broadPeak.sorted");
  foreach my $file (@files) {
    my $output = "test".$opc;
    my $cmd = "perl /tools/trick/intersectFiles.pl -o $input -m $file -oiend 2 -miend 2 -count >$output";
    system($cmd);
    $input = $output;
    $opc += 1;
  }
}
