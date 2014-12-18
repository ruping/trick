use strict;

my $gene_list = shift;

my $database = shift;

open IN, "$gene_list";
my %genes;
while ( <IN> ){
   chomp;
   $genes{$_} = '';
}
close IN;
my $np = scalar(keys %genes);

open DA, "$database";
while ( <DA> ){
  next if /^pathway\t/;
  chomp;
  my ($name, $id, $source, $symbols) = split /\t/;
  my @symbols = split (',', $symbols);
  my $white;
  my $nw = scalar(@symbols);
  my $pw = 0;
  foreach my $s (@symbols){
    if (exists ($genes{$s})){
       $white .= $s.',';
       $pw += 1;
    }
  }
  print "$_\t$white\t$nw\t$np\t$pw\n";
}
close DA;
