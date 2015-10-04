use strict;

my $samheader = shift;
my $chrom = shift;

if ($chrom !~ /^[cC][hH][rR]/){
  $chrom = 'chr'.$chrom;
}

open IN, "$samheader";
while ( <IN> ){
  if ( $_ =~ /^\@SQ\tSN\:(\S+)\tLN\:\d+/ ) {
    my $chr = $1;
    if ($chr !~ /^[cC][hH][rR]/){
      $chr = 'chr'.$chr;
    }
    if ( $chrom eq $chr ) {
      print "$_";
    }
  } else {
    print "$_";
  }
}
close IN;
