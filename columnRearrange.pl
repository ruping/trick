use strict;
use Data::Dumper;

my $file = shift;
my $order = shift;

my @order = split(/\,/, $order);

open IN, "$file";
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  if ($order ne ''){
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
          if ($cols[$i] =~ /^TCGA-.+?$/ or $cols[$i] =~ /^AC\d+/){
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
          push(@order, $i+1) if ($cols[$i] =~ /^TCGA/ or $cols[$i] =~ /^AC\d+/);
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
