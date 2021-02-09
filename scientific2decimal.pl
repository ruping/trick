use strict;

my $file = shift;
my $column = shift;
my @columns = split (',', $column);

open IN, "$file";
while ( <IN> ){
  chomp;
  next if /^[\#]?[Cc][Hh][Rr]\t/;
  my @cols = split /\t/;
  foreach my $column (@columns){
    $cols[$column] = sprintf("%.10g", $cols[$column]);
  }
  printf("%s\n", join("\t", @cols));
}
close IN;
