use strict;

open IN, shift;
my %data;
while ( <IN> ){
  chomp;
  if (/^[\#]?chr/){
    print "$_\n";
  }
  else {
    my @cols = split /\t/;
    $data{$cols[0]}{$cols[1]} = $_;
  }
}
close IN;

foreach my $chr (sort {$a cmp $b} keys %data){
  foreach my $pos (sort {$a <=> $b} keys %{$data{$chr}}){
     print "$data{$chr}{$pos}\n";
  }
}
