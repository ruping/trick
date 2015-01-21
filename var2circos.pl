use strict;
use Data::Dumper;
use Getopt::Long;

my $file;
my $expr;
my $prefix;
my $sample;
my $rnaseq;

GetOptions (
           "file|f=s"       => \$file,            #variant result filename
           "expr|e=s"       => \$expr,            #comma seperated indexes
           "sample|s=s"     => \$sample,
           "rnaseq|r=s"     => \$rnaseq,
           "help|h"         => sub{
                               print "usage: $0 generate circos files from variants and expression\n\nOptions:\n\t--file\t\tthe filename of variants\n";
                               print "\t--expr\t\texpression table of rpkm\n";
                               print "\t--rnaseq\t\trnaseq data id\n";
                               print "\t--sample\t\tsample in request\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );


my %expr;
open EX, "$expr";
while ( <EX> ){
  next if /^#/;
  chomp;
  my ($ens, $count, $RPKM, $genetype, $geneName) = split /\t/;
  $expr{$geneName} = sprintf("%.3f", log2($RPKM + 0.001));
}
close EX;


#chr     pos     id      ref     alt     TB002   TB002maf        TB002d  TB003   TB003maf        TB003d  function.       somatic germline        rep     sc      bp      TB001   TB001d  cmeanav cmedianav       clin
my @name;
my %colnames;

open IN, "$file";
open OUT1, ">$file\.mutGene\.label";
open OUT2, ">$file\.mutDepthTotal";
open OUT3, ">$file\.mutDepthMut";
open OUT4, ">$file\.mutDepthTotalRNA";
open OUT5, ">$file\.mutDepthMutRNA";
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  if ($_ =~ /^[\#]?chr\t/){
    $_ =~ s/^#//;
    @name = @cols;
    for (my $i = 0; $i <= $#name; $i++) {
      $colnames{$name[$i]} = $i;
    }
    next;
  } else {
    #for (my $i == 0; $i <= $#cols; $i++) {
    #  if ($cols[$i])
    #}
    $cols[$colnames{'chr'}] =~ s/^chr//;
    $cols[$colnames{'chr'}] = 'hs'.$cols[$colnames{'chr'}];
    if ($cols[$colnames{'function.'}] =~ /;geneName=(.+?);/){
      my $geneName = $1;
      print OUT1 "$cols[$colnames{'chr'}] $cols[$colnames{'pos'}] $cols[$colnames{'pos'}] $geneName\n";
    }
    my $depthTotal = sprintf("%.3f", log2($cols[$colnames{$sample.'d'}] + 1));
    print OUT2 "$cols[$colnames{'chr'}] $cols[$colnames{'pos'}] $cols[$colnames{'pos'}] $depthTotal\n";
    my $depthMut = sprintf("%.3f", log2($cols[$colnames{$sample.'d'}] + 1)*$cols[$colnames{$sample.'maf'}]);
    print OUT3 "$cols[$colnames{'chr'}] $cols[$colnames{'pos'}] $cols[$colnames{'pos'}] $depthMut\n";
    my $depthTotalRNA = sprintf("%.3f", log2($cols[$colnames{$rnaseq.'d'}] + 1));
    print OUT4 "$cols[$colnames{'chr'}] $cols[$colnames{'pos'}] $cols[$colnames{'pos'}] $depthTotalRNA\n";
    my $depthMutRNA = sprintf("%.3f", log2($cols[$colnames{$rnaseq.'d'}] + 1)*$cols[$colnames{$rnaseq}]);
    print OUT5 "$cols[$colnames{'chr'}] $cols[$colnames{'pos'}] $cols[$colnames{'pos'}] $depthMutRNA\n";
  }
}
close IN;
close OUT1, "";
close OUT2, "";
close OUT3, "";
close OUT4, "";
close OUT5, "";

sub log2 {
  my $n = shift;
  return log($n)/log(2);
}

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}
