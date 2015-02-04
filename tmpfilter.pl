use strict;
use Data::Dumper;

my $file = shift;
my $list = shift;

open IN, "$list";
my %vcf;
while ( <IN> ){
  chomp;
  my $name;
  if ($_ =~ /(AC\d+)[^a-zA-Z0-9]/){
    $name = $1;
  }
  $vcf{$name} = $_;
}
close IN;

print STDERR Dumper(\%vcf);
