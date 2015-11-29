use strict;

my $header = shift;
my $sample = shift;

open IN, "$header";
while ( <IN> ){
  if ($_ =~ /\@RG/){
    $_ =~ s/\tSM\:[A-Z0-9a-z\.\_\-]+/\tSM\:$sample/;
  }
  print "$_";
}
close IN;
