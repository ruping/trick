use strict;

my $file = shift;
my $dummy = "SRP";

open IN, "$file";
while ( <IN> ){
  chomp;
  s/[\n\s]$//;
  my @cols = split /\t/;
  push(@cols, $dummy);
  printf("%s\n",join("\t",@cols));
}
close IN;

exit 22;
