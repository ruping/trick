use strict;

my $nameMapping = shift;
my $file = shift;
my $prefix = shift;
my $columns = shift;

my @columns = split(',', $columns);
my %needc;
foreach my $col (@columns){
  $needc{$col} = '';
}

open MAP, "$nameMapping";
my %mapping;
while ( <MAP> ){
  chomp;
  my @cols = split /\t/;
  $mapping{$cols[0]} = $cols[2];
}
close MAP;

open IN, "$file";
while ( <IN> ) {
  chomp;
  my @cols = split /\t/;
  if (/^[\#]?chr\t/){  #header
    for (my $i = 0; $i <= $#cols; $i++){
      if (exists ($mapping{$cols[$i]})){
        my $newName = $mapping{$cols[$i]};
        $cols[$i] = $newName;
      }
    }
  }

  else{ #normal lines
    for (my $i = 0; $i <= $#cols; $i++) {

      next if (!exists($needc{$cols[$i]}));

      my @names = split(',', $cols[$i]);
      foreach my $name (@names){
        $name =~ /^($prefix\d+)/;
        my $rname = $1;
        my $newName = $mapping{$rname};
        $name =~ s/^$prefix\d+/$newName/;
      }
      $cols[$i] = join(',', @names);
    }
  }
  printf("%s\n", join("\t", @cols));
}
close IN;
