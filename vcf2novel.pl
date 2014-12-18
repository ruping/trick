use strict;

my $file = shift;
my $whichNovel = shift;
my $mode = "novel";
if ($whichNovel =~ /^(.+?)\,(.+?)$/){
  $mode = "common";
}
my $indexNovel = 1;

print STDERR "mode: $mode\n";

open IN, "$file";
while ( <IN> ) {
  chomp;
  if (/^#/) {
    if (/^#CHROM\t/) {
      my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $GT1, $GT2) = split /\t/;
      if ($mode eq "novel") {
        if ($GT2 eq $whichNovel) {
          $indexNovel = 2;
        }
      }
    }
    print "$_\n";
    next;
  }

  my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $GT1, $GT2) = split /\t/;

  $GT1 =~ /^(\d\/\d)\:([\.0-9]+\:)?\d+\,(\d+)\:(\d+)\:/;
  my $genot1 = $1;
  my $ad1 = $2;
  my $dp1 = $3;
  $GT2 =~ /^(\d\/\d)\:([\.0-9]+\:)?\d+\,(\d+)\:(\d+)\:/;
  my $genot2 = $1;
  my $ad2 = $2;
  my $dp2 = $3;

  if ($mode eq "novel") {
    if ($indexNovel == 1) {
      if ($genot1 eq '0/1' and $genot2 eq '0/0' and $ad2 == 0 and $dp2 >= 8) {
        print "$_\n";
      }
    } elsif ($indexNovel == 2) {
      if ($genot2 eq '0/1' and $genot1 eq '0/0' and $ad1 == 0 and $dp1 >= 8) {
        print "$_\n";
      }
    }
  }

  if ($mode eq "common") {
    if ($genot1 eq $genot2) {
       print "$_\n";
    }
  }

}
close IN;
