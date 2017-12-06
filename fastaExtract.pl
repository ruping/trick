use strict;

my $genome = shift;
my $chrs = shift;

my @chrs = split(',', $chrs);
my %chrs;
foreach my $chr (@chrs) {
  $chrs{$chr} = '';
}

open HS, "$genome";
my %genome = ();
my $chr = undef;
while ( <HS> ) {

  if (/^>(\w+).*?\n$/) {
    $chr = $1;
  }
  if ( exists($chrs{$chr}) ) {
    print "$_";
  }

}
close HS;

exit 0;
