use strict;

my $fastq = shift;

open IN, "$fastq";

while ( <IN> ){

  $_ =~ /^@(\S+)\s+([12])\S+$/;
  my $identifier = $1;
  my $add = $2;
  print ">$identifier\/$add\n";

  $_ = <IN>;
  print "$_";

  $_ = <IN>;
  $_ = <IN>;

}

close;
