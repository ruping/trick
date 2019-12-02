use Data::Dumper;
use Parallel::ForkManager;

@chrs = qw(1 2 3 4 5 6 7 8 9 10);

my $chrBatches = partitionArray(\@chrs, 3);

foreach my $chrBatch (@{$chrBatches}) {
  my $manager = Parallel::ForkManager->new($options{'threads'});
  my $processedChroms = "chromosome ";
  foreach my $chrom (@{$chrBatch}) {
    $manager->start and next;
    my $cmd = "echo $chrom >$chrom";
    system( $cmd );
    $manager->finish;
    $processedChroms .= $chrom.',';
  }
  $manager->wait_all_children;
  print STDERR "$processedChroms have been processed!\n";
}

sub partitionArray {
  my ($arr, $N) = @_;

  my @res;
  my $i = 0;

  while ($i + $N-1 <= $#$arr) {
    push @res, [@$arr[$i .. $i+$N-1]];
    $i += $N;
  }

  if ($i <= $#$arr) {
    push @res, [@$arr[$i .. $#$arr]];
  }
  return \@res;
}
