use strict;
use File::Glob ':glob';
use Data::Dumper;

my $dir = shift;

my @tissue = qw(Rectal_Mucosa Rectal_Smooth_Muscle Duodenum_Mucosa Lung Pancreatic_Islets Colonic_Mucosa);

foreach my $tissue (@tissue) {
  my @files = bsd_glob("$dir/$tissue/*.narrowPeak.sorted");
  foreach my $file (@files) {
    my $out = $file;
    $out =~ s/narrowPeak.sorted/narrowPeak.sorted.sig/;
    open OUT, ">$out";
    open IN, "$file";
    while ( <IN> ) {
      chomp;
      next if /^#/;
      my ($chr, $start, $end, $name, $score, $strand, $signal, $pValue, $qValue) = split /\t/;
      next if ($qValue <= 1.3); #qvalue threshold 0.05
      printf OUT ("%s\n", join("\t", $chr, $start, $end));
    }
    close IN;
    close OUT;
  }
}




