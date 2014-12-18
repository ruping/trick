open HS, shift;

my %genome = ();

my $chr = "SRP";

while ( <HS> ){

   if (/^>/) {
     close OUT unless ($chr eq 'SRP');

     $_ =~ /^(\S+)$/;
     $chr = $1;

     #unless ($chr eq 'mitochondrion') {
     #  $chr = 'chr'.$chr;
     #}

     #if ($chr eq 'mitochondrion'){
     #  $chr = 'chrM';
     #}

     open OUT, ">$chr\.fa";
     print OUT "$_";
   } else {
     print OUT "$_";
   }

 }

close HS;
