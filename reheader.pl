use strict;

my $file = shift;

my $add = shift;
my @add = split(",", $add);
my $realadd = join("\t", @add);

open IN, "$file";
while ( <IN> ){
  chomp;
  if (/^#/){
    print "$_\t$realadd\n";
  } else {
    print "$_\n";
  }
}
close IN;
