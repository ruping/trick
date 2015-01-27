use strict;

my $files = shift;

my @files = split(/\,/, $files);

#########################if paola table exist########################################
my $paola = shift;
my %paola;
my @paolaNames;
if ($paola ne '') {
   open IN, "$paola";
   while ( <IN> ) {
    chomp;
    next if ($_ =~ /^[\@\#]/);
    my @cols = split /\t/;
    if ($_ =~ /^gene\t/) {
      @paolaNames = @cols;
      next;
    } else {
      my $gene;
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($paolaNames[$i] eq 'gene') {
          $gene = $cols[$i];
          next;
        }
        $paola{$gene}{$paolaNames[$i]} = $cols[$i];
      }                           #iterator
    }                             #else
  }
  close IN;
}
############################paola table###############################################


my %type2int;
$type2int{'snv'} = 2;
$type2int{'snv&raregermline'} = 3;
$type2int{'indel'} = 4;
$type2int{'indel&raregermline'} = 5;
$type2int{'fusion'} = 6;
$type2int{'cnva'} = 1;
$type2int{'cnvd'} = -1;

my @ileum = qw(AC444 AC445 AC446 AC516 AC517 AC518 AC519);
my %ileum;
foreach my $il (@ileum) {
  $ileum{$il} = '';
}
my @rectum = qw(AC3 AC439 AC440 AC441 AC442 AC443 AC447 AC525 AC526 AC527 AC528 AC529 AC530 AC531 AC532 AC533 AC546 AC548 AC580 AC637 AC653 AC668);
my %rectum;
foreach my $rec (@rectum){
  $rectum{$rec} = '';
}

my %result;
foreach my $file (@files) {
  my $type;
  if ($file =~ /snv/){
    $type = 'snv';
    if ($file =~ /germline/){
      $type .= '&raregermline';
    }
  } elsif ($file =~ /indel/){
    $type = 'indel';
    if ($file =~ /germline/){
      $type .= '&raregermline';
    }
  } elsif ($file =~ /fusion/){
    $type = 'fusion';
  } elsif ($file =~ /copynumber/){
    $type = 'cnv';
  }

  open IN, "$file";
  my @name;
  while ( <IN> ) {
    chomp;
    next if ($_ =~ /^[\@\#]/);
    my @cols = split /\t/;
    if ($_ =~ /^gene\t/) {
      @name = @cols;
      next;
    } else {
      my $gene;
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($name[$i] eq 'gene'){
          $gene = $cols[$i];
        }
        if ($name[$i] =~ /(AC\d+)/) {
          my $sample = $1;
          if ($type eq 'cnv') {
            if ($cols[$i] eq 'NA') {
              #do nothing
            } elsif ($cols[$i] > 1) {
              $type = $type.'a';
            } elsif ($cols[$i] <= 1) {
              $type = $type.'d';
            }
          }
          $result{$gene}{$sample}{$type} = $cols[$i];
        }                         #it is sample
      }                           #iterator
    }                             #else
  }
  close IN;
}

if ($paola ne ''){
  goto PAOLA;
}

print "gene";
foreach my $sample (@rectum){
  print "\t$sample";
}
foreach my $sample (@ileum){
  print "\t$sample";
}
print "\trectum\tileum\n";

foreach my $gene (keys %result) {
  print "$gene";
  my $rectum = 0;
  my $ileum = 0;
  foreach my $sample (@rectum) {     #rectum
    my $changed = 0;
    my $vars;
    foreach my $type (%{$result{$gene}{$sample}}) {
      if ($result{$gene}{$sample}{$type} > 0 and $result{$gene}{$sample}{$type} ne 'NA') {
         $changed = 1;
         #$vars .= $type.","
         $vars .= $type2int{$type};
      }
    }
    if ($vars eq ''){
      $vars = 0;
    }
    $rectum += $changed;
    print "\t$vars";
  }
  foreach my $sample (@ileum) {     #ileum
    my $changed = 0;
    my $vars;
    foreach my $type (sort {my $tia = $type2int{$a}; my $tib = $type2int{$b}; $a <=> $b} keys %{$result{$gene}{$sample}}) {
      if ($result{$gene}{$sample}{$type} > 0){
         $changed = 1;
         #$vars .= $type.","
         $vars .= $type2int{$type};
      }
    }
    if ($vars eq ''){
      $vars = 0;
    }
    $ileum += $changed;
    print "\t$vars";
  }
  print "\t$rectum\t$ileum\n";
}

exit 0;

PAOLA:

print "gene";
foreach my $paolaName (@paolaNames) {
  print "\t$paolaName" unless ($paolaName eq 'gene');
}
print "\n";

foreach my $gene (sort {$a cmp $b} keys %paola) {
  next if $gene =~ /\"/;

  $paola{$gene}{'sumPz'} = '';
  $paola{$gene}{'SOMMUT'} = '';

  my %tmptype;
  foreach my $paolaName (@paolaNames) {
    if ($paolaName =~ /AC\d+/) {
      if ($paola{$gene}{$paolaName} =~ /(-?1)/) {my $dt = $1; $tmptype{$dt}++;}
      if ($result{$gene}{$paolaName} ne '') {
        $paola{$gene}{$paolaName} =~ s/^0$//;
        $paola{$gene}{$paolaName} =~ s/[23]//g;
        foreach my $type (sort {my $tia = $type2int{$a}; my $tib = $type2int{$b}; $tia <=> $tib} keys %{$result{$gene}{$paolaName}}) {
           $paola{$gene}{$paolaName} .= $type2int{$type} if $result{$gene}{$paolaName}{$type} > 0;
           $tmptype{$type2int{$type}}++ if $result{$gene}{$paolaName}{$type} > 0;
        }
        $paola{$gene}{$paolaName} = 0 if ($paola{$gene}{$paolaName} eq '');
      } else {   #for old stuff
        while ($paola{$gene}{$paolaName} =~ /([23])/g) {
           my $dt = $1;
           $tmptype{$dt}++;
        }
      }
      $paola{$gene}{'sumPz'} += ($paola{$gene}{$paolaName} == 0)? 0:1;
    } #each sample

    if ($paolaName eq 'SOMMUT') {
      foreach my $tmpt (sort {$a <=> $b} keys %tmptype) {
        $paola{$gene}{'SOMMUT'} .= $tmpt;
      }
      if ($paola{$gene}{'SOMMUT'} eq ''){
        $paola{$gene}{'SOMMUT'} = 0;
      }
    }

    if ($paolaName eq 'gene'){
      print "$gene";
    } else {
      print "\t$paola{$gene}{$paolaName}";
    }
  } #paolaName [col names]
  print "\n";
}

