use strict;

my $sv = shift;
#my $sns = shift;
#my @sns = split(",", $sns);

my %colindex;
my %colnames;
open IN, "$sv";
while ( <IN> ) {
    chomp;
    my @cols = split /\t/;
    if ($_ =~ /^AnnotSV/ ) {  #header
	for (my $i = 0; $i <= $#cols; $i++) {
            $colindex{$cols[$i]} = $i;
            $colnames{$i} = $cols[$i];
	}
        #printf("%s\n", join("\t", "Hugo_Symbol","Chromosome","Start_Position","End_position","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode"));
    } else {

        next if ($cols[$colindex{'location'}] eq 'txStart-txEnd');
        next if ($cols[$colindex{'AnnotSV type'}] eq 'full');
        $cols[$colindex{'SV type'}] =~ s/^DEL$/Deletion/;
        $cols[$colindex{'SV type'}] =~ s/^DUP$/Duplication/;
        $cols[$colindex{'SV type'}] =~ s/^BND$/Translocation/;
        $cols[$colindex{'SV type'}] =~ s/^INV$/Inversion/;

        my $startidx = $colindex{'FORMAT'}+1;
        my $endidx = $colindex{'AnnotSV type'}-1;
        for my $idx ($startidx..$endidx) {
          my $sample = $colnames{$idx};;
          if ($cols[$colindex{$sample}] =~ /^[01]\/1/) {  #called
              printf("%s\n", join("\t", $cols[$colindex{'Gene name'}], $cols[$colindex{'SV chrom'}], $cols[$colindex{'SV start'}], $cols[$colindex{'SV end'}], $cols[$colindex{'SV type'}],
                                        'SV', $cols[$colindex{'REF'}], $cols[$colindex{'ALT'}], $sample));
          }

        }
    }
}
close IN;
