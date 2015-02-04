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


open IN, "$file";
my %data;
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  $data{$cols[2]}{$cols[0]}{$cols[1]} = '';
}
close IN;

foreach my $sample (keys %data){
  my $calling = $vcf{$sample};
  open OUT, ">./tmp";
  foreach my $chr (sort keys %{$data{$sample}}){
    foreach my $pos (sort {$a <=> $b} keys %{$data{$sample}{$chr}}){
      print OUT "$chr\t$pos\t$sample\n";
    }
  }
  close OUT;
  my $cmd = "perl /ifs/home/c2b2/ac_lab/rs3412/tools/trick/intersectFiles.pl -o ./tmp -m $calling >>./result";
  RunCommand($cmd);
}

sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}
