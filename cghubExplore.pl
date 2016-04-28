use strict;

open IN, shift;
my %need;
while ( <IN> ) {
  chomp;
  my $id;
  my $sampleId;
  my $analysisId;
  my $study;
  my $bamfile;
  my $baifile;
  if ($_ =~ /^\s+(Analysis\s+\d+)$/) {
    $id = $1;
  } else {
    next if $id !~ /^Analysis/;
    if ($_ =~ /^\s+legacy_sample_id\s+\:\s+(.+?)$/) {
      $sampleId = $1;
      next if $sampleId =~ /^CCLE/;
    } elsif ($_ =~ /^\s+analysis_id\s+\:\s+(.+?)$/) {
      $analysisId = $1;
    } elsif ($_ =~ /^\s+study\s+\:\s+(.+?)$/) {
      $study = $1;
    } elsif ($_ =~ /^\s+filename\s+\:\s+(.+?\.bam)$/) {
      $bamfile = $1;
      next if $bamfile =~ /^CCLE/;
    } elsif ($_ =~ /^\s+filename\s+\:\s+(.+?\.bai)$/) {
      $baifile = $1;
    } elsif ($_ =~ /^\s+analysis\_data\_uri/) { #save
      print "$sampleId\n";
      $need{$sampleId}{'aid'} = $analysisId;
      $need{$sampleId}{'bam'} = $bamfile;
      $need{$sampleId}{'bai'} = $baifile;
    }
  }
}
close IN;

foreach my $sampleId (sort keys %need) {
  my $bam = $need{$sampleId}{'bam'};
  my $aid = $need{$sampleId}{'aid'};
  print "$aid\t$bam\n";
}

exit 0;
