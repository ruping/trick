use strict;

my $original = shift;
my $realmaf = shift;

open OR, "$original";
my %OR;
while ( <OR> ){
   chomp;
   next if /^#/;
   #chr    pos     id      ref     alt     AC3     AC439   AC440   AC441   AC442   AC443   AC444   AC445   AC446   AC447   AC516   AC517   AC518   AC519   AC525   AC526   AC527   AC528   AC529   AC530   AC531   AC532   AC533   AC546   AC548   AC580   function
   my @cols = split /\t/;
   my $coor = $cols[0].':'.$cols[1];
   $OR{$coor} = $cols[$#cols];
}
close OR;
print STDERR "$original loaded\n";

open RM, "$realmaf";
while ( <RM> ){
  if ($_ =~ /^#/){
     next;
  }
  chomp;
  my ($chr,$pos,$id,$ref,$alt,$AC1,$AC1d,$AC3,$AC3d,$AC439,$AC439d,$AC440,$AC440d,$AC441,$AC441d,$AC442,$AC442d,$AC443,$AC443d,$AC444,$AC444d,$AC445,$AC445d,$AC446,$AC446d,$AC447,$AC447d,$AC516,$AC516d,$AC517,$AC517d,$AC518,$AC518d,$AC519,$AC519d,$AC525,$AC525d,$AC526,$AC526d,$AC527,$AC527d,$AC528,$AC528d,$AC529,$AC529d,$AC530,$AC530d,$AC531,$AC531d,$AC532,$AC532d,$AC533,$AC533d,$AC546,$AC546d,$AC548,$AC548d,$AC580,$AC580d) = split /\t/;
  my $coor = $chr.':'.$pos;
  my $rectum = 0;
  my $ileum;

  #if (($AC3>0 and $AC439>0 and $AC440>0 and $AC441>0 and $AC443>0 and $AC447>0 and $AC525>0 and $AC526>0 and $AC527>0 and $AC528>0 and $AC529>0 and $AC530>0 and $AC531>0 and $AC548>0) and ($AC444 == 0 and $AC445==0 and $AC446==0 and $AC516==0 and $AC517==0 and $AC518==0 and $AC519==0)){
    # print "$_\t$OR{$coor}\n";
  #}

  $rectum += 1 if ($AC3>0);
  $rectum += 1 if ($AC439>0);
  $rectum += 1 if ($AC440>0);
  $rectum += 1 if ($AC441>0);
  $rectum += 1 if ($AC443>0);
  $rectum += 1 if ($AC447>0);
  $rectum += 1 if ($AC525>0);
  $rectum += 1 if ($AC526>0);
  $rectum += 1 if ($AC527>0);
  $rectum += 1 if ($AC528>0);
  $rectum += 1 if ($AC529>0);
  $rectum += 1 if ($AC530>0);
  $rectum += 1 if ($AC531>0);
  $rectum += 1 if ($AC548>0);

  $ileum += 1 if ($AC444>0);
  $ileum += 1 if ($AC445>0);
  $ileum += 1 if ($AC446>0);
  $ileum += 1 if ($AC516>0);
  $ileum += 1 if ($AC517>0);
  $ileum += 1 if ($AC518>0);
  $ileum += 1 if ($AC519>0);

  if ($rectum >= 12 and $ileum <= 1){
    print "$_\t$OR{$coor}\n";
  }

}
close RM;
