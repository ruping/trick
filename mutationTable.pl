use strict;

my $mutation = shift;
my $task = shift;

my @type1 = qw(AC2maf AC3maf AC4maf AC53maf AC54maf AC55maf AC56maf AC57maf AC565maf AC566maf AC567maf) if ($task eq 'arj');
my %type1;
foreach my $sample (@type1){
  $type1{$sample} = 0;
}
my @type2 = qw(AC582maf AC583maf) if ($task eq 'arj');
my %type2;
foreach my $sample (@type2){
  $type2{$sample} = 0;
}

open IN, "$mutation";
my @colnames;
while ( <IN> ) {
  chomp;
  if ($_ =~ /^chr/) {
    @colnames = split /\t/;
    for(my $c = 0; $c <= $#colnames; $c++){
      if ($c == 0){
        print "$colnames[$c]";
      } elsif ($colnames[$c] eq 'id'){
        print "\tlink\tid";
      } else {
        print "\t$colnames[$c]";
      }
    }
    print "\tcloneType\n" if ($task eq 'arj');
    print "\n" if ($task eq 'cohort');
  } else {
    my @cols = split /\t/;
    my $chr;
    my $pos;
    my $ucsc;
    my $nt1 = 0;
    my $nt2 = 0;
    my $nt1nz = 0;
    my $nt2nz = 0;
    my $primary = 0;
    my $cloneType = "somatic";
    for (my $i = 0; $i <= $#colnames; $i++) {
      if ($colnames[$i] eq 'chr') {
        $chr = $cols[$i];
        print "$chr";
      } elsif ($colnames[$i] eq 'pos') {
        $pos = $cols[$i];
        $ucsc = "\=HYPERLINK\(\"http\:\/\/genome\.ucsc\.edu\/cgi\-bin/hgTracks\?db\=hg19\&position\=chr$chr\%3A$pos\-$pos\", \"UCSC\"\)";
        print "\t$pos\t$ucsc";
      } elsif ($colnames[$i] eq 'id' or $colnames[$i] eq 'ref' or $colnames[$i] eq 'alt') {
        print "\t$cols[$i]";
      } elsif ( $colnames[$i] =~ /AC\d+maf/ ) {
        print "\t$cols[$i]";
        if (exists ($type1{$colnames[$i]})){
          $nt1++ if ($cols[$i] >= 0.1);
          $nt1nz++ if ($cols[$i] > 0);
        } elsif (exists ($type2{$colnames[$i]})){
          $nt2++ if ($cols[$i] >= 0.1);
          $nt2nz++ if ($cols[$i] > 0);
        }
        if ( $colnames[$i] eq 'AC565maf' and $cols[$i] == 0 ){
          $primary += 0.5;
        }
      } elsif ( $colnames[$i] =~ /^AC\d+d$/ ) {
        if ( $colnames[$i] eq 'AC565d' and $cols[$i] >= 5 ){
          $primary += 0.5;
        }
        print "\t$cols[$i]";
      } else {
        print "\t$cols[$i]";
      }
    }
    if (($nt1 >= 1 and $nt1nz >= 5) and ($nt2 > 0 and $nt2nz > 1)){
      $cloneType .= '.both'
    } elsif (($nt1 <= 1 and $nt1nz < 3) and ($nt2 > 0 and $nt2nz > 1)){
      $cloneType .= '.type2';
    } elsif (($nt2 == 0 and $nt2nz == 0) and $nt1 >= 2){
      $cloneType .= '.type1';
    }
    if ($primary == 1){
      $cloneType .= '.metastasis';
    }
    print "\t$cloneType\n" if ($task eq 'arj');
    print "\n" if ($task eq 'cohort');
  }
}
close IN;
