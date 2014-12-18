use strict;

my $file = shift;
my $type = shift;

open IN, "$file";
while ( <IN> ){
  chomp;
  if ($_ =~ /^#/){
    print "$_\n";
    next;
  }
  my ($chr, $pos, $id, $ref, $alt, $qual, $pass, $info, $format, $sample, $blood) = split /\t/;
  if ($info =~ /INDEL/ || $info =~ /Indel/){
     print "$_\n" if ($type eq 'indel');
  } else {
     print "$_\n" if ($type eq 'snv');
  }
}
close IN;
