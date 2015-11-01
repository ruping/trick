use strict;

my $gene_list = shift;

my $database = shift;

open IN, "$gene_list";
my %genes;
while ( <IN> ) {
  chomp;
  next if /^[\#]?gene\t/;
  my @cols = split /\t/;
  $genes{$cols[0]} = '';
}
close IN;
my $np = scalar(keys %genes);

open DA, "$database";
my %pathgenes;
while ( <DA> ) {
  next if /^pathway\t/;
  chomp;
  my ($name, $id, $source, $symbols) = split /\t/;
  my @symbols = split (',', $symbols);
  my $white;
  my $nw = scalar(@symbols);
  my $pw = 0;
  foreach my $s (@symbols){
    $pathgenes{$s} = '';
    if (exists ($genes{$s})){
       $white .= $s.',';
       $pw += 1;
    }
  }
  print "$_\t$white\t$nw\t$np\t$pw\n";
}
close DA;

my $npathgenes = scalar(keys %pathgenes);
print STDERR "$npathgenes\n";
