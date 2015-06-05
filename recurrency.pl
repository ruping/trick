use strict;
use Data::Dumper;

my $file = shift;
my $maf = shift;
my $somaticInfo = shift; #if for somatic judgement

my %somatic;
my %germline;  #may have multiple tumors
if ($somaticInfo ne '' and -s "$somaticInfo") {

  open IN, "$somaticInfo";
  while ( <IN> ){
    chomp;
    s/[\s\n]$//;
    my @columns = split /\t/;
    my $tumor = $columns[0];
    my $normal = $columns[1];

    $somatic{$tumor} = $normal;
    push(@{$germline{$normal}}, $tumor) if $normal ne 'undef';
  }
  close IN;
  print STDERR Dumper (\%somatic);
  print STDERR Dumper (\%germline);
}


my @rectum = qw(AC57maf AC439maf AC440maf AC441maf AC443maf AC447maf AC525maf AC526maf AC527maf AC528maf AC529maf AC530maf AC531maf AC532maf AC533maf AC546maf AC548maf AC580maf AC637maf AC653maf AC668maf);
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
while ( <IN> ) {
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
    } elsif ($maf =~ /trace/) {
      print "$_\ttrace\n";
    } elsif ($maf eq 'founds') {
      print "$_\tfounds\n";
    } elsif ($maf eq 'somatic') {
      print "$_\tsomatic\tgermline\n";
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
    } elsif ($maf == 1) {    #get real maf
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
    } elsif ($maf =~ /trace/) {  #trace sample
      my $trace = '';
      foreach my $sample (@all) {
        my $samp = $sample;
        if ($samp ne 'AC3maf') {
          $samp =~ s/maf$//;
          if ($cols[$colnames{$samp}] > 0){
              $trace .= "$samp,"
          }
        } else {  #AC3maf
          if ($cols[$colnames{$samp}] >= 0.1) {
            my $vard = sprintf("%.1f", $cols[$colnames{$samp}]*$cols[$colnames{$samp}+1]);
            if ($vard >= 2) {
              $samp =~ s/maf$//;
              $trace .= "$samp,";
            }
          }
        }
      }
      $trace =~ s/,$//;
      if ($trace ne '') {
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
        print "$_\t$ftrace\n" if ($maf eq 'trace');
        print "$_\t$trace\n" if ($maf eq 'traceall');
      }
    } elsif ($maf eq 'founds') {  #trace all samples
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
    } elsif ($maf eq 'somatic'){  #find somatic ones
      my %tumor;
      my %blood;
      my %nonblood;
      my %unknown;
      for ( my $i = 0; $i <= $#cols; $i++ ) {
        if ($colindex{$i} =~ /^(.+?)maf$/) {
          my $samp = $1;
          my $maf = $cols[$i];
          my $depth = $cols[$i+1];
          my $vard = sprintf("%.1f", $maf*$depth);

          if ($cols[$colnames{'coor'}] == 55545257) {
            print STDERR "$samp\t";
            print STDERR "$maf\t";
            print STDERR "$depth\t";
            print STDERR "$vard\n";
          }

          if (exists $somatic{$samp}) { #it is tumor
             $tumor{$samp} = $maf if ($vard >= 2 and $maf >= 0.1);
          } elsif (exists $germline{$samp}) { #it is blood
            foreach my $ct (@{$germline{$samp}}) {
              if ( $maf == 0 and $depth >= 10 ) {
                $nonblood{$ct} = '';
              } elsif ( $maf == 0 and $depth < 10 ) {
                $unknown{$ct} = '';
              } elsif ($vard >= 2) {
                $blood{$ct} = $maf;
              }
            }
          } #it is blood
        } #maf
      } #each column

      if ($cols[$colnames{'coor'}] == 55545257){
        print STDERR Dumper (\%tumor);
        print STDERR Dumper (\%blood);
        print STDERR Dumper (\%nonblood);
        print STDERR Dumper (\%unknown);
      }

      my $soma = 'NA';
      my $germ = 'NA';
      foreach my $tumorSamp (keys %tumor) {
        my $stype = 'NA';
        if (exists $nonblood{$tumorSamp}) {
          $stype = 'good';
          $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';
        } elsif (exists($blood{$tumorSamp})) {
          if ($blood{$tumorSamp} <= 0.1 and $tumor{$tumorSamp}/$blood{$tumorSamp} >= 4) {
            $stype = 'doubt';
            $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';
          } elsif (exists($unknown{$tumorSamp})) {
            $stype = 'undef';
            $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';
          } else {
            $germ = ($germ eq 'NA')? $tumorSamp.',':$germ.$tumorSamp.',';
          }
        } else { #unknown ones, all germ
          $germ = ($germ eq 'NA')? $tumorSamp.',':$germ.$tumorSamp.',';
        }
      }
      print "$_\t$soma\t$germ\n";
    } #somatic
  }
}
close IN;
