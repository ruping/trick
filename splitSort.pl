use strict;
use File::Glob ':glob';
use Data::Dumper;

#broadPeak
#chrom - Name of the chromosome (or contig, scaffold, etc.).
#chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
#chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
#name - Name given to a region (preferably unique). Use '.' if no name is assigned.
#score - Indicates how dark the peak will be displayed in the browser (0-1000).
#strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
#signalValue - Measurement of overall (usually, average) enrichment for the region.
#pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
#qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned. 


#narrowPeak
#chrom - Name of the chromosome (or contig, scaffold, etc.).
#chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
#chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
#name - Name given to a region (preferably unique). Use '.' if no name is assigned.
#score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
#strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
#signalValue - Measurement of overall (usually, average) enrichment for the region.
#pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
#qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
#peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called. 


my $dir = shift;
my $suffix = shift;
my $instead = shift;
my $type = shift;

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
  printf OUT ("%s\n", join("\t", '#chr', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pValue', 'qValue')) if ($type eq 'broadPeak');
  printf OUT ("%s\n", join("\t", '#chr', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pValue', 'qValue', "peak")) if ($type eq 'narrowPeak');
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
