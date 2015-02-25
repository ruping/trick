use strict;

my $file = shift;
my $maf = shift;


my @rectum = qw(AC3maf AC439maf AC440maf AC441maf AC443maf AC447maf AC525maf AC526maf AC527maf AC528maf AC529maf AC530maf AC531maf AC532maf AC533maf AC546maf AC548maf AC580maf AC637maf AC653maf AC668maf);
my @ileum = qw(AC444maf AC445maf AC446maf AC516maf AC517maf AC518maf AC519maf);
my @primary = qw(AC532maf AC533maf AC546maf AC580maf AC668maf);
my @blood = qw(AC1maf AC547maf AC581maf AC669maf);
my @all = ();
push(@all, @rectum);
push(@all, @ileum);
push(@all, "AC442maf");


open IN, "$file";
my %colnames;
my %colindex;
while ( <IN> ){
  chomp;
  if (/^[\#]?chr\t/){
    #it is header
    my @cols = split /\t/;
    for(my $i = 0; $i <= $#cols; $i++) {
      $colnames{$cols[$i]} = $i;
      $colindex{$i} = $cols[$i];
    }
    if ($maf eq '') {
      print "$_\tfounds\tfounds.rectum\tfounds.ileum\tfounds.primary\n";
    } elsif ($maf == 1) {
      print "$_\tmaf\n";
    } elsif ($maf eq 'trace') {
      print "$_\ttrace\n";
    } elsif ($maf eq 'founds') {
      print "$_\tfounds\n";
    }
  } else {
    my @cols = split /\t/;
    if ($maf eq ''){
      my $founds = 0;
      my $foundsRectum = 0;
      my $foundsIleum = 0;
      my $foundsPrimary = 0;
      foreach my $rec ( @rectum ) {
        if ($cols[$colnames{$rec}] >= 0.1) {
          my $vard = sprintf("%.1f", $cols[$colnames{$rec}]*$cols[$colnames{$rec}+1]);
          if ($vard >= 2) {
            $foundsRectum++;
            $founds++;
          }
        }
      }

      foreach my $ile ( @ileum ) {
        if ($cols[$colnames{$ile}] >= 0.1) {
          my $vard = sprintf("%.1f", $cols[$colnames{$ile}]*$cols[$colnames{$ile}+1]);
          if ($vard >= 2) {
            $foundsIleum++;
            $founds++;
          }
        }
      }

      foreach my $pri (@primary) {
        if ($cols[$colnames{$pri}] >= 0.1) {
          my $vard = sprintf("%.1f", $cols[$colnames{$pri}]*$cols[$colnames{$pri}+1]);
          if ($vard >= 2) {
            $foundsPrimary++;
          }
        }
      }

      print "$_\t$founds\t$foundsRectum\t$foundsIleum\t$foundsPrimary\n";
    } elsif ($maf == 1) {
      my $maf = 0;
      my $sampleCounts = scalar(@all);
      foreach my $sample (@all) {
        if ($cols[$colnames{$sample}] >= 0.1) {
          my $vard = sprintf("%.1f", $cols[$colnames{$sample}]*$cols[$colnames{$sample}+1]);
          if ($vard >= 2) {
            $maf += $cols[$colnames{$sample}];
          }
        }
      }
      $maf = sprintf("%.6f",$maf/$sampleCounts);
      print "$_\t$maf\n";
    } elsif ($maf eq 'trace') {
      my $trace = '';
      foreach my $sample (@all) {
        $sample =~ s/maf$// unless ($sample eq 'AC3maf');
        if ($sample ne 'AC3maf') {
          if ($cols[$colnames{$sample}] > 0){
              $trace .= "$sample,"
          }
        } else {
          if ($cols[$colnames{$sample}] >= 0.1) {
            my $vard = sprintf("%.1f", $cols[$colnames{$sample}]*$cols[$colnames{$sample}+1]);
            if ($vard >= 2) {
              $sample =~ s/maf$//;
              $trace .= "$sample,";
            }
          }
        }
      }
      $trace =~ s/,$//;
      if ($trace ne ''){
        my @trace = split(/\,/, $trace);
        my $ftrace = '';
        if (scalar(@trace) > 1) {
          if ($trace[0] eq 'AC3') {
            $ftrace = $trace[1];
          } else {
            $ftrace = $trace[0];
          }
        } else {
          $ftrace = $trace[0];
        }
        print "$_\t$ftrace\n";
      }
    } elsif ($maf eq 'founds') {
      my $founds = 0;
      for (my $i = 0; $i <= $#cols; $i++){
        if ($colindex{$i} =~ /maf$/) {
          if ($cols[$i] >= 0.1){
            my $vard = sprintf("%.1f", $cols[$i]*$cols[$i+1]);
            if ($vard >= 2) {
              $founds++;
            }
          }
        } #maf
      } #each column
      print "$_\t$founds\n";
    } #get founds
  }
}
close IN;
