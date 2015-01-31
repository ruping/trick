use strict;

my $file = shift;

my @rectum = qw(AC3maf AC439maf AC440maf AC441maf AC443maf AC447maf AC525maf AC526maf AC527maf AC528maf AC529maf AC530maf AC531maf AC532maf AC533maf AC546maf AC548maf AC580maf AC637maf AC653maf AC668maf);
my @ileum = qw(AC444maf AC445maf AC446maf AC516maf AC517maf AC518maf AC519maf);
my @primary = qw(AC532maf AC533maf AC546maf AC580maf AC668maf);
my @blood = qw(AC1maf AC547maf AC581maf AC669maf);

open IN, "$file";
my %colnames;
while ( <IN> ){
  chomp;
  if (/^[\#]?chr\t/){
    #it is header
    my @cols = split /\t/;
    for(my $i = 0; $i <= $#cols; $i++){
      $colnames{$cols[$i]} = $i;
    }
    print "$_\tfounds\tfounds.rectum\tfounds.ileum\tfounds.primary\n";
  } else {
    my @cols = split /\t/;
    my $founds = 0;
    my $foundsRectum = 0;
    my $foundsIleum = 0;
    my $foundsPrimary = 0;
    foreach my $rec (@rectum){
      if ($cols[$colnames{$rec}] >= 0.1){
        my $vard = sprintf("%.1f", $cols[$colnames{$rec}]*$cols[$colnames{$rec}+1]);
        if ($vard >= 2){
          $foundsRectum++;
          $founds++;
        }
      }
    }

    foreach my $ile (@ileum){
      if ($cols[$colnames{$ile}] >= 0.1){
        my $vard = sprintf("%.1f", $cols[$colnames{$ile}]*$cols[$colnames{$ile}+1]);
        if ($vard >= 2){
          $foundsIleum++;
          $founds++;
        }
      }
    }

    foreach my $pri (@primary){
      if ($cols[$colnames{$pri}] >= 0.1){
        my $vard = sprintf("%.1f", $cols[$colnames{$pri}]*$cols[$colnames{$pri}+1]);
        if ($vard >= 2){
          $foundsPrimary++;
        }
      }
    }

    print "$_\t$founds\t$foundsRectum\t$foundsIleum\t$foundsPrimary\n";
  }
}
close IN;
