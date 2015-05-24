use strict;
use Data::Dumper;
use Getopt::Long;

my $file;
my $order;
my $prefix;

GetOptions (
           "file|f=s"       => \$file,             #filename
           "order|o=s"      => \$order,            #comma seperated indexes
           "prefix|p=s"     => \$prefix,
           "help|h"         => sub{
                               print "usage: $0 rearrange columns according to your order or for maf rearrange\n\nOptions:\n\t--file\t\tthe filename to be reordered\n";
                               print "\t--order\t\tcomma seperated indexes\n";
                               print "\t--prefix\tthe prefix of samples' names\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );


my @order = split(/\,/, $order);

open IN, "$file";
while ( <IN> ) {
  chomp;
  my @cols = split /\t/;
  if ($order ne ''){ #order
    my $print;
    foreach my $index (@order) {
      $print .= $cols[$index]."\t";
    }
    $print =~ s/\t$//;
    print "$print\n";
  } else { #order not specified, for variants maf table
    if (/^#/) {
      my %samples;
      my $secondstart = 0;

      for (my $i = 0; $i <= $#cols; $i++){
        if ($cols[$i] eq 'bp') {
          $secondstart = $i+1;
          push(@order, $i);
          last;
        } else {
          push(@order, $i);
          if ($cols[$i] =~ /^TCGA-.+?$/ or $cols[$i] =~ /^$prefix\d+/){
            $samples{$cols[$i]} = $i;
          }
        } #else
      } #for

      my %inserted;
      for (my $i = $secondstart; $i <= $#cols; $i++) {
        if ( exists($samples{$cols[$i]}) ) {
          my $insertpos = $samples{$cols[$i]};
          my $offset;
          foreach my $sample ( sort { my $sa = $samples{$a}; my $sb = $samples{$b}; $sa <=> $sb } keys %samples) {
             my $rank = $samples{$sample};
             last if ($insertpos == $rank);
             if ( exists($inserted{$sample}) ) {
                $offset += 2;
             }
          }
          $insertpos += $offset+1;
          splice(@order, $insertpos, 0, $i, $i+1);
          $inserted{$cols[$i]} = '';
          $cols[$i] .= 'maf';
        } elsif ($cols[$i] !~ /d$/) {
          push(@order, $i);
          push(@order, $i+1) if ($cols[$i] =~ /^TCGA/ or $cols[$i] =~ /^$prefix\d+/);
        }
      }

    } #if it is the header

    my $print;
    foreach my $index (@order){
      $print .= $cols[$index]."\t";
    }
    $print =~ s/\t$//;
    print "$print\n";

  } #column not specified
}
close IN;
