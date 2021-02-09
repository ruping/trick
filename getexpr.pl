use strict;
use File::Glob ':glob';
use File::Basename;

open NEED, shift;
my %need;
while (<NEED>){
  chomp;
  $need{$_} = '';
}
close NEED;

#open STILL, ">>/ifs/scratch/c2b2/ac_lab/rs3412/no1/net/still_need2.txt";
my $bdir = "/ifs/scratch/c2b2/ac_lab/rs3412/no1/brca/ruping/rnaseq/bams/";
my $outdir = "/ifs/scratch/c2b2/ac_lab/rs3412/no1/brca/ruping/rnaseq/expr/";

#my @allbams = bsd_glob("$bdir/*.bam");
#foreach my $allbam (@allbams){
#   my $tmp = basename($allbam);
#   $tmp =~ s/\.bam$//;
#   if (! exists $need{$tmp}){
#     print "$tmp\n";
#   }
#}

foreach my $need (sort keys %need) {

  if (-e "$bdir/$need\.bam") { #do counting

   my $bam = "$bdir/$need\.bam";
   unless (-e "$outdir/$need\_gencode.expr") {
      my $cmd = "/ifs/home/c2b2/ac_lab/rs3412/tools/TRUP/bin/Rseq_bam_reads2expr --region /ifs/data/c2b2/ac_lab/rs3412/annotation/hg19/hg19.gencode/gencode.v14.annotation.gene.bed12 --mapping $bam >$outdir/$need\_gencode.expr";
      unless (-e "$outdir/$need\_gencode.expr"){
        #system($cmd);
      }
   }

   unless (-e "$outdir/$need\_gencode.expr.sorted"){
      my $cmd = "sort -k 1,1d -k 2,2n $outdir/$need\_gencode.expr >$outdir/$need\_gencode.expr.sorted";
      system($cmd);
   }

   my $entrez_file = "$outdir/$need\_gencode.expr.sorted.entrez";
   my $gencode_file = "$outdir/$need\_gencode.expr.sorted";
   unless (-e $entrez_file) {
    my $cmd = "perl entrez_bed.pl /ifs/data/c2b2/ac_lab/rs3412/annotation/hg19/hg19.biomart.txt $gencode_file >$entrez_file";
    print STDERR "$cmd\n";
    system($cmd);
   }

  } else {
    print STDERR "$bdir/$need\.bam not found\n";
  }

}
