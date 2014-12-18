use strict;

my $fastq = shift;
my $add = shift;

open IN, "$fastq";

while ( <IN> ){

  chomp;
  $_ =~ /@(\S+)\s+(\S+)\s+(\S+)$/;
  my $identifier = $2;
  print "\@$identifier\/$add\n";

  $_ = <IN>;
  print "$_";

  $_ = <IN>;
  $_ =~ /@(\S+)\s+(\S+)\s+(\S+)$/;
  $identifier = $2;
  print "\+$identifier\/$add\n";

  $_ = <IN>;
  print "$_";

}

close;
