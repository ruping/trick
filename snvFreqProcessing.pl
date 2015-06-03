use strict;

open IN, shift;
my %colindex;
while ( <IN> ) {
  chomp;
  my @cols = split /\t/;
  if ($_ =~ /^[\#]?chr\t/){  #header rows
    for (my $i = 0; $i <= $#cols; $i++){
     $colindex{$cols[$i]} = $i;
    }
  } else {  #data rows
    my $id = $cols[$colindex{'id'}];
    my $esp = $cols[$colindex{'ESP6500SI.snv.vcf.INFO.MAF'}];
    my $dbsnp = $cols[$colindex{'dbsnp142.snv.vcf.INFO.CAF'}];
    my $somatic = $cols[$colindex{'somatic'}];
    my $founds = $cols[$colindex{'founds'}];

    
  }
}
close IN;
