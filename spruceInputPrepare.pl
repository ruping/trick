use strict;
use List::Util qw[max];

open IN, shift;
my $print;
my $m;
my $n;
my @allsampindex;
my @allmutindex;
while ( <IN> ){
  chomp;
  if (/^sample\_index/){   #header
    my $header = $_;
    $header =~ s/\tseg$//;
    $header = '#'.$header;
    $print .= $header."\n";
  } else {
    my ($sample_index, $sample_label, $character_label, $character_index, $vaf_lb, $vaf_mean, $vaf_ub, $x1, $y1, $mu1, $x2, $y2, $mu2, $seg) = split /\t/;
    my $buffer = join("\t",$sample_index, $sample_label, $character_label, $character_index, $vaf_lb, $vaf_mean, $vaf_ub, $x1, $y1, $mu1, $x2, $y2, $mu2);
    $buffer =~ s/\tNA//g;
    $print .= $buffer."\n";
    push(@allsampindex, $sample_index);
    push(@allmutindex, $character_label);
  }
}
close IN;

$m = max(@allsampindex)+1;
$n = max(@allmutindex)+1;
print "$m #m\n$n #n\n$print";
