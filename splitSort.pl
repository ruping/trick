use strict;
use File::Glob ':glob';
use Data::Dumper;

my $dir = shift;
my $suffix = shift;
my $instead = shift;

my @files = bsd_glob("$dir/*.$suffix");
print Dumper(\@files);
my $totalN = scalar(@files);
print "$totalN\n";

foreach my $input (@files){
  open IN, "$input";
  my %data;
  while (<IN>) {
    chomp;
    my @cols = split /\t/;
    $data{$cols[0]}{$cols[1]} = $_;
    #next if /^track name/;
    #my ($ID, $chr, $start, $end, $gene, $rank, $super, $H3K27ac, $control) = split /\,/;
    #$data{$chr}{$start}{'end'} = $end;
    #$data{$chr}{$start}{'id'} = $ID;
    #$data{$chr}{$start}{'gene'} = $gene;
    #$data{$chr}{$start}{'rank'} = $rank;
    #$data{$chr}{$start}{'super'} = $super;
    #$data{$chr}{$start}{'H3K27ac'} = $H3K27ac;
    #$data{$chr}{$start}{'control'} = $control;
  }
  close IN;

  (my $output = $input) =~ s/$suffix$/$instead/;
  open OUT, ">$output";
  #printf OUT ("%s\n", join("\t", '#chr', 'start', 'end', 'id', 'gene', 'rank', 'super', 'H3K27ac', 'control'));
  foreach my $chr (sort keys %data) {
    foreach my $start (sort {$a <=> $b} keys %{$data{$chr}}) {
      #my $end = $data{$chr}{$start}{'end'};
      #my $id = $data{$chr}{$start}{'id'};
      #my $gene = $data{$chr}{$start}{'gene'};
      #my $rank = $data{$chr}{$start}{'rank'};
      #my $super = $data{$chr}{$start}{'super'};
      #my $H3K27ac = $data{$chr}{$start}{'H3K27ac'};
      #my $control = $data{$chr}{$start}{'control'};
      #printf OUT ("%s\n", join("\t", $chr, $start, $end, $id, $gene, $rank, $super, $H3K27ac, $control));
      print OUT "$data{$chr}{$start}\n";
    }
  }
  close OUT;
}
