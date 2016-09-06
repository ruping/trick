use strict;

my $ms = shift;

my @mappingStats = split (',', $ms);

my @stats;
my %stats;
foreach my $mappingStats (@mappingStats){
  open IN, "$mappingStats";
  while ( <IN> ) {
    chomp;
    $_ =~ /^(\w+\:)(\s+)(\d+)$/;
    my $statn = $1;
    push (@stats, $statn);
    my $space = $2;
    my $count = $3;
    $stats{$statn}{'space'} = $space;
    $stats{$statn}{'count'} += $count;
  }
  close IN;
}

foreach my $statn (@stats){
  printf("%s\n",join("", $statn, $stats{$statn}{'space'}, $stats{$statn}{'count'}));
}

exit 0;
