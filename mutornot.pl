use strict;


my @name;
my %found;

#my $input = $ARGV[0];
#open IN, "$input";
while ( <> ) {
  chomp;
  my @cols = split /\t/;
  if ($_ =~ /^chr\t/) {
    @name = @cols;
    next;
  } else {
    for (my $i = 0; $i <= $#cols; $i++) {
      if ($name[$i] =~ /(AC\d+)maf/) {
        my $name = $1;
        if ($cols[$i] >= 0.1) {
          $found{$name}++;      #found
        }
      }                         #it is maf
    }                           #iterator
  }                             #else
}                               #while
#close IN;

for (my $i = 0; $i <= $#name; $i++) {
  if ($name[$i] =~ /(AC\d+)maf/) {
    my $name = $1;
    if (exists($found{$name})){
       print "$name\t$found{$name}\n";
    } else {
       print "$name\t0\n";
    }
  }                             #it is maf
}                               #iterator
