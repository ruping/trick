use strict;
use Data::Dumper;

my $wig = shift;
my $region = shift;
my $ratio = shift;

my @region = split(',', $region);
my %region;
foreach my $r (@region) {
  $r =~ /^(.+?)\:(\d+)\-(\d+)$/;
  my $chr = $1;
  my $start = $2;
  my $end = $3;
  $region{$chr} = [$start, $end];
}
print STDERR Dumper(\%region);

open W, "$wig";
my ($chr, $start, $step, $span, $cpos);

while ( <W> ) {
  chomp;
  if (/^fixedStep\schrom\=(.+?)?\sstart\=(\d+)\sstep\=(\d+)\sspan\=(\d+)/) {
    $chr = $1;
    $start = $2;
    $step = $3;
    $span = $4;
    $cpos = $start;
    print "$_\n";
  } else {
    my $count = $_;
    $cpos += $step;
    if (exists($region{$chr})) {
      my $start = $region{$chr}->[0];
      my $end = $region{$chr}->[1];
      if ($cpos >= $start and $cpos <= $end) { #modi!
        $count = round($count*($ratio));
        print STDERR "modi: $chr\t$cpos\t$count\n";
      }
    }
    print "$count\n";
  }
}
close W;

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}
