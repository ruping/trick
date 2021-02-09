use strict;
use Data::Dumper;

my $wigs = shift;
my @wigs = split(',', $wigs);

my @blocks;
my %merged;

my $wigNum = 1;
foreach my $wig (@wigs) {
  open W, "$wig";
  my ($chr, $start, $step, $span, $cpos);
  my $block = '';
  my $line = 0;
  while ( <W> ) {
    chomp;
    if (/^fixedStep\schrom\=(.+?)?\sstart\=(\d+)\sstep\=(\d+)\sspan\=(\d+)/) {
      $chr = $1;
      $start = $2;
      $step = $3;
      $span = $4;
      $cpos = $start;
      $line = 1;
      $block = "$_";
      push(@blocks, $block) if ($wigNum == 1);
    } else {
      my $count = $_;
      $merged{$block}{$line} += $count;
      $line++;
    }
  }
  close W;
  $wigNum++;
}

for (my $i = 0; $i <= $#blocks; $i++) {
  print "$blocks[$i]\n";
  foreach my $line (sort {$a <=> $b} keys %{$merged{$blocks[$i]}}) {
    print "$merged{$blocks[$i]}{$line}\n";
  }
}

exit 0;
