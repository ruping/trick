use strict;
use Data::Dumper;

my $need = shift;

my $dir = shift;

my @need = split(",", $need);


my %rna2dna;
my @arj = qw(AC10 AC13 AC16 AC35 AC38 AC41 AC47 AC50 AC501 AC502 AC584 AC585 AC586);
my %arj;
foreach my $arj (@arj) {
  $arj{$arj} = '';
}

open IN, "/ifs/scratch/c2b2/ac_lab/rs3412/no1/net/dna2rna.mapping";
while ( <IN> ) {
  chomp;
  next if /^dna\t/;
  my ($dna, $rna, $avail) = split /\t/;
  if ($avail eq 'Y') {
    if (exists ($arj{$rna})){
      $rna2dna{$rna} = 'AC3';
    } else {
      $rna2dna{$rna} = $dna;
    }
  }
}
close IN;

#print Dumper (\%rna2dna);

my %fusions;

foreach my $sample (@need) {

  (my $sampleName = $sample) =~ s/\_001//;
  my $sampleDNA = $rna2dna{$sampleName};
  my $ff = $dir."/".$sample."/05_FUSION/$sampleName\.fusion\.report";

  if ($sampleName eq 'AC252' or $sampleName eq 'AC263' or $sampleName eq 'AC281' or $sampleName eq 'AC283') {
    $ff = $dir."../30M_SE/".$sample."/05_FUSION/$sampleName\.fusion\.report";
  }

  if (! -e $ff){
    print STDERR "$ff does not exist!!!\n";
    exit 22;
  }

  open FF, "$ff" || die "could not open $ff\n";
  while ( <FF> ) {
    chomp;
    my ($fusion1, $fusion2, $bp1, $bp2, $dir, $rep, $sc, $type, $strand, $cov1, $cov2, $cov3, $cov4, $cov5, $transcripts) = split /\t/;

    my $gene1 = '';
    my $gene2 = '';
    if ($fusion1 =~ /^\S+\((.+?)\)$/){
      $gene1 = $1;
    } else {
      $gene1 = $fusion1;
    }

    if ($fusion2 =~ /\S+\((.+?)\)/){
      $gene2 = $1;
    } else {
      $gene2 = $fusion2;
    }

    #for LRIG1#########################################
    if ($gene1 eq 'FRYL' and $gene2 eq '') {
      if ($bp2 =~ /(chr\w+)\:(\d+)/){
        if ($1 eq 'chr3' and $2 > 66612400 and $2 < 66612500) {
           $gene2 = 'LRIG1';
        }
      }
    }
    if ($gene2 eq 'FRYL' and $gene1 eq '') {
      if ($bp1 =~ /(chr\w+)\:(\d+)/){
        if ($1 eq 'chr3' and $2 > 66612400 and $2 < 66612500) {
           $gene1 = 'LRIG1';
        }
      }
    }
    ###################################################

    next if ($rep =~ /R/ and $sc =~ /C/);                     #skip both rep and sc
    next if ($fusion1 =~ /KCNMB2/ and $fusion2 =~ /KCNMB2/);  #skip KCNMB2
    next if ($fusion1 =~ /TMPRSS3/ or $fusion2 =~ /TMPRSS3/); #skip TMPRSS3
    next if ($cov5 < 1 or $cov3 < 2);                         #skip low spanning reads
    next if ($cov5 < 2 and $sc eq 'CC');                      #skip low spanning reads
    next if ($cov4/$cov5 >= 10 and $cov5 <= 3);               #high duplication
    next if ($cov4/$cov5 >= 20 and $cov5 < 10);               #high duplication
    next if ((($gene2 ne 'IGR' and $gene1 =~ /^$gene2/) or ($gene1 ne 'IGR' and $gene2 =~ /^$gene1/)) and $sc eq 'CC');  #remaining family
    next if (($gene1 =~ /^IGK[JV]/ or $gene2 =~ /^IGK[JV]/) and $sc eq 'CC');                                            #IGG
    next if (($gene1 eq 'IGR' and $gene2 =~ /^RP\d+\-\d+/) or ($gene2 eq 'IGR' and $gene1 =~ /^RP\d+\-\d+/));            #lnRNA IGR junk

    if ($gene1 ne ""){
      $fusions{$gene1}{$sampleDNA} = $cov5;
    }
    if ($gene2 ne ""){
      $fusions{$gene2}{$sampleDNA} = $cov5;
    }
  }

  close FF;

}

my %dna;
foreach my $rna (keys %rna2dna){
  $dna{$rna2dna{$rna}} = '';
}
print "gene";
foreach my $sample ( sort keys %dna ){
   print "\t$sample";
}
print "\n";

foreach my $gene ( sort {$a cmp $b} keys %fusions ) {
  print "$gene";
  foreach my $sampleName ( sort keys %dna ) {
    if ( $fusions{$gene}{$sampleName} ne '' ) {
       print "\t$fusions{$gene}{$sampleName}";
    } else {
       print "\t0";
    }
  }
  print "\n";
}
