use strict;

open IN, shift;
my $oldChr = 'SRP';
my $oldPos = 0;
my $oldCall = "";
while ( <IN> ){
  chomp;
  if ($_ =~ /^chr\t/){
    print "$_\n";
    next;
  }
  my @cols = split /\t/;
  my $chr = $cols[0];
  my $pos = $cols[1];

  if ($chr ne $oldChr) {
    if ($oldChr ne 'SRP'){
      print "$oldCall\n";
    }
    $oldChr = $chr;
    $oldPos = $pos;
    $oldCall = $_;
  } else {
    my $dis = $pos - $oldPos;
    if ($dis > 10){
      if ($oldPos != 0){
        print "$oldCall\n";
      }
      $oldPos = $pos;
      $oldCall = $_;
    } else {
      #do nothing
    }
  }

}
close IN;
