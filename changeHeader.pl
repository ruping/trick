use strict;

my $file = shift;
my $change = shift;

my @changes = split(',', $change);
my %changes;
for(my $i = 0; $i <= $#changes; $i = $i+2) {
  $changes{$changes[$i]} = $changes[$i]+1;
}

open IN, $file;
while ( <IN> ) {
  chomp;
  if (/^[#]?chr\t/){
    my @cols = split /\t/;
    foreach my $c (keys %changes){
      $cols[$c] = $changes{$c};
    }
    printf("%s\n", join("\t", @cols));
  } else {
    print "$_\n";
  }
}
close IN;
