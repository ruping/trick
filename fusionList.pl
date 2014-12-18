use strict;

my %rna2dna;
my @arj = qw(AC10 AC11 AC12 AC13 AC14 AC15 AC16 AC17 AC18 AC19 AC20 AC21 AC35 AC36 AC37 AC38 AC39 AC40 AC41 AC42 AC43 AC47 AC48 AC49 AC50 AC51 AC52 AC501 AC502 AC584 AC585 AC586);
my %arj;
foreach my $arj (@arj) {
  $arj{$arj} = '';
}

my @cohort = qw(AC241 AC242 AC243 AC245 AC246 AC261 AC274 AC294 AC373 AC252 AC263 AC281 AC283 AC452 AC455 AC534 AC535 AC538 AC577);
my %cohort;
foreach my $cohort (@cohort) {
  $cohort{$cohort} = '';
}

open IN, "/cygdrive/f/dna2rna.mapping";
my %rna2dna;
while ( <IN> ) {
  chomp;
  next if /^dna\t/;
  my ($dna, $rna, $avail, $sampleID) = split /\t/;
  if ($avail eq 'Y') {
     $rna2dna{$rna} = $dna."\t".$sampleID;
  }
}
close IN;

my $ff = shift;

printf("%s\n", join("\t", qw(group RNA_id DNA_id sample_id seqDepth gene1 gene2 breakpoint1 breakpoint2 direction repetitive selfChain type strand Span_Reads_nR Encompass_Reads_nR Total_Reads_nR Span_R Span_Score Assembled_contig_id)));

open FF, "$ff" || die "could not open $ff\n";
my %fusions;
while ( <FF> ) {
  chomp;
  my ($column1, $fusion2, $bp1, $bp2, $dir, $rep, $sc, $type, $strand, $cov1, $cov2, $cov3, $cov4, $cov5, $transcripts) = split /\t/;

  next if ($rep =~ /R/ and $sc =~ /C/); #skip both rep and sc
  $column1 =~ /^(\d+M\_[PS]E)\//;
  my $seqDepth = $1;
  $column1 =~ /05\_FUSION\/(AC\d+)\.fusion\.report\:(.+?)$/;
  my $rna = $1;
  my $dnaandsample = $rna2dna{$rna};
  if ($dnaandsample eq ''){
    $dnaandsample = "\t";
  }
  my $identifier = $rna;
  my $fusion1 = $2;

  next if ($fusion1 =~ /KCNMB2/ and $fusion2 =~ /KCNMB2/); #skip KCNMB2
  next if ($fusion1 =~ /TMPRSS3/ or $fusion2 =~ /TMPRSS3/); #skip TMPRSS3
  next if ($cov5 < 1);          #skip low spanning reads

  my $gene1 = '';
  my $gene2 = '';
  if ($fusion1 =~ /^\S+\((.+?)\)$/) {
    $gene1 = $1;
  } elsif ($fusion1 eq 'IGR') {
    $gene1 = 'IGR';
  }

  if ($fusion2 =~ /\S+\((.+?)\)/) {
    $gene2 = $1;
  } elsif ($fusion2 eq 'IGR') {
    $gene2 = 'IGR';
  }

  $bp1 =~ /(chr\w+)\:(\d+)/;
  my $chr1 = $1;
  my $coor1 = $2;
  $bp2 =~ /(chr\w+)\:(\d+)/;
  my $chr2 = $1;
  my $coor2 = $2;
  my $dis = abs($coor1-$coor2);
  #print STDERR "$dis\n";

  if (($chr1 eq $chr2) and $dis < 1000000){
     $type = "Read_Th";
  }


  #for LRIG1#########################################
  if ($gene1 eq 'FRYL' and $gene2 eq 'IGR') {
     if ($chr2 eq 'chr3' and $coor2 > 66612400 and $coor2 < 66612500) {
        $gene2 = 'LRIG1';
     }
  }
  if ($gene2 eq 'FRYL' and $gene1 eq 'IGR') {
     if ($chr1 eq 'chr3' and $coor1 > 66612400 and $coor1 < 66612500) {
       $gene1 = 'LRIG1';
     }
  }
  ###################################################

  my $value = sprintf("%s\n", join("\t",($rna,$dnaandsample,$seqDepth,$gene1,$gene2,$bp1,$bp2,$dir,$rep,$sc,$type,$strand,$cov1,$cov2,$cov3,$cov4,$cov5,$transcripts)));
  push(@{$fusions{$identifier}}, $value);
}

close FF;


foreach my $rnasample (@arj){
  foreach my $value (@{$fusions{$rnasample}}){
    print "patient0\t$value";
  }
}

foreach my $rnasample (@cohort){
  foreach my $value (@{$fusions{$rnasample}}){
    print "net_cohort\t$value";
  }
}

foreach my $rnasample (sort keys %fusions) {
  next if exists($arj{$rnasample});
  next if exists($cohort{$rnasample});
  foreach my $value (@{$fusions{$rnasample}}){
    print "other\t$value";
  }
}
