use strict;

my $id = "SRP";
my $sampleId = "SRP";
my $analysisId = "SRP";
my $study = "SRP";
my $bamfile = "SRP";
my $baifile = "SRP";
my $center = "SRP";
my $plf = "SRP";
my $lib = "SRP";
my %need;

open IN, shift;

while ( <IN> ) {
  chomp;
  if ($_ =~ /^\s+(Analysis\s+\d+)$/) {
    $id = $1;
  } else {
    next if $id !~ /^Analysis/;
    if ($_ =~ /^\s+legacy_sample_id\s+\:\s+(.+?)$/) {
      $sampleId = $1;
    } elsif ($_ =~ /^\s+analysis_id\s+\:\s+(.+?)$/) {
      $analysisId = $1;
    } elsif ($_ =~ /^\s+study\s+\:\s+(.+?)$/) {
      $study = $1;
    } elsif ($_ =~ /^\s+center_name\s+\:\s+(\S+)$/) {
      $center = $1;
    } elsif ($_ =~ /^\s+filename\s+\:\s+(.+?\.bam)$/) {
      $bamfile = $1;
    } elsif ($_ =~ /^\s+filename\s+\:\s+(.+?\.bai)$/) {
      $baifile = $1;
    } elsif ($_ =~ /^\s+platform\s+\:\s+(\S+)$/){
      $plf = $1;
    } elsif ($_ =~ /^\s+library_strategy\s+\:\s+(\S+)$/){
      $lib = $1;
    } elsif ($_ =~ /^\s+analysis\_data\_uri/) { #save
      next if ($sampleId =~ /^CCLE/ | $bamfile =~ /^CCLE/);
      $need{$sampleId}{'aid'} = $analysisId;
      $need{$sampleId}{'bam'} = $bamfile;
      $need{$sampleId}{'bai'} = $baifile;
      $need{$sampleId}{'cen'} = $center;
      $need{$sampleId}{'plf'} = $plf;
      $need{$sampleId}{'lib'} = $lib;
    }
  }
}
close IN;

foreach my $sampleId (sort keys %need) {
  my $ind;
  if ($sampleId =~ /^(TCGA\-[A-Z0-9]+\-[A-Z0-9]+)/){
    $ind = $1;
  }
  my $bam = $need{$sampleId}{'bam'};
  my $aid = $need{$sampleId}{'aid'};
  my $cen = $need{$sampleId}{'cen'};
  my $plf = $need{$sampleId}{'plf'};
  my $lib = $need{$sampleId}{'lib'};
  print "$ind\t$sampleId\t$bam\t$cen\t$plf\t$lib\n";
}

exit 0;
