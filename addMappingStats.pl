use strict;

my $ms = shift;

my @mappingStats = split (',', $ms);

my @stats;
my %stats;
my $fileCount = 1;
foreach my $mappingStats (@mappingStats){
  open IN, "$mappingStats";
  while ( <IN> ) {
    chomp;
    $_ =~ /^(\w+\:)(\s+)(\d+)$/;
    my $statn = $1;
    push (@stats, $statn) if ($fileCount == 1);
    my $space = $2;
    my $count = $3;
    $stats{$statn}{'space'} = $space;
    $stats{$statn}{'count'} += $count;
  }
  close IN;
  $fileCount ++;
}

foreach my $statn (@stats) {
  printf("%s\n",join("", $statn, $stats{$statn}{'space'}, $stats{$statn}{'count'}));
}

exit 0;
