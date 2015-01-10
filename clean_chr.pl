my @chrs = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM);
my %chrs;

foreach my $id (@chrs){
  $chrs{$id} = '';
}

my $file = shift;
my $remove = shift;
#print "#$file\n";

open IN, "$file";
while ( <IN> ){

   if (/^#/){
     print "$_";
     next;
   }

   if ($remove ne ''){
     $_ =~ s/^chr//;
     print "$_";
     next;
   }

   $_ =~ /^(\S+)/;
   my $chr_now = $1;
   my $chr_now_full = 'chr'.$chr_now;

   if ( $chr_now eq 'chrMT' ) {
      $_ =~ s/^chrMT/chrM/;
      print "$_";
      next;
   }
   elsif ($chr_now_full eq 'chrMT'){
      $_ =~ s/^MT/M/;
      print "chr$_";
      next;
   }


   if (exists $chrs{$chr_now_full}) {
      print "chr$_";
   }
   elsif (exists $chrs{$chr_now}){
      print "$_";
   }

}
close IN;
