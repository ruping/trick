use strict;

my $old_SC = shift;
my $chain_gff= shift;

#bin    score   tName   tSize   tStart  tEnd    qName   qSize   qStrand qStart  qEnd    id
# 0       1       2       3       4      5        6       7       8       9      10     11

open IN, "$chain_gff";
my %chain;
while ( <IN> ){
  next if /^#/;
  chomp;
  my @cols = split /\t/;
  $cols[8] =~ /L=(\d+)$/;
  my $linec = $1;
  my $hash = {'chr'=>$cols[0], 'start'=>$cols[3], 'end'=>$cols[4]};
  push (@{$chain{$linec}}, $hash);
}
close IN;

my $line = 1;
open IN, "$old_SC";
while ( <IN> ){
  next if /^#/;
  chomp;
  my @cols = split /\t/;
  if (exists($chain{$line}) and scalar(@{$chain{$line}}) == 2){ #do the chain change
    $cols[2] = ${$chain{$line}}[0]->{'chr'};
    $cols[4] = ${$chain{$line}}[0]->{'start'};
    $cols[5] = ${$chain{$line}}[0]->{'end'};
    $cols[6] = ${$chain{$line}}[1]->{'chr'};
    $cols[9] = ${$chain{$line}}[1]->{'start'};
    $cols[10] = ${$chain{$line}}[1]->{'end'};
    printf("%s\n", join("\t", @cols));
  }
  $line++;
}
close IN;
