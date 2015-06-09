use strict;

my $bin = shift;
my $transfac = shift;
my $fasta = shift;

my %res;

open IN, "$transfac";
undef $/;                    #redefine seperator
$/ = "//";
my $tf;
my @tf;
my $count = 0;
while ( <IN> ) {
  $_ =~ s/\/\///;
  $_ =~ s/^\n//;
  if ($_ =~ /^VV.+?\nXX\n$/) {
    next;
  }

  if ($_ =~ /\nID\s+(\S+)\n/) {
    $tf = $1;
    push (@tf, $tf);
  }
  if ($count <= 1000) {
    open OUT, ">transfacTmp";
    print OUT "$_";
    close OUT;
    my $cmd = "$bin/trap transfacTmp $fasta >resultTmp";
    RunCommand($cmd, 0, 0);

    open RES, "<resultTmp";
    undef $/;
    $/ = "\n";
    while ( <RES> ){
      chomp;
      my @cols = split /\t/;
      if ($cols[1] >= 1) {
        $res{$cols[0]}{$tf} = $cols[1];
      }
    }
    close RES;
    undef $/;
    $/ = "//";

    $count++;
  } else {
    last;
  }
}
close IN;

my $size = scalar(keys %res);
print STDERR "size: $size\n";

#print "#coor";
#foreach my $tf (@tf) {
#    print "\t$tf";
#}
#print "\n";

foreach my $region (keys %res) {
  print "$region";
  foreach my $tf (sort keys %{$res{$region}}) {
    if ($res{$region}{$tf} =~ /\d/){
      print "\t$tf\=$res{$region}{$tf}";
    }
    #elsif (! defined $res{$region}{$tf}) {
    #  print "\t0";
    #}
  }
  print "\n";
}


exit;


sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}
