use strict;

my $file = shift;
my $type = shift;
my $cohort = shift;

my @ileum = qw(AC444 AC445 AC446 AC516 AC517 AC518 AC519);
my %ileum;
foreach my $il (@ileum) {
  $ileum{$il} = '';
}
my @rectum = qw(AC3 AC439 AC440 AC441 AC442 AC443 AC447 AC525 AC526 AC527 AC528 AC529 AC530 AC531 AC532 AC533 AC546 AC548 AC580 AC637 AC653 AC668);
my %rectum;
foreach my $rec (@rectum) {
  $rectum{$rec} = '';
}

open IN, "$file";
my @name;
my %result;
my @samples;

while ( <IN> ) {
  chomp;
  next if ($_ =~ /^[\@\#]/);
  my @cols = split /\t/;
  if ($_ =~ /^chr\t/) {
    @name = @cols;
    next;
  } else {
    my @genes = &grepGene($_);
    for (my $i = 0; $i <= $#cols; $i++) {
      if ($name[$i] =~ /(AC\d+)maf/) {
        my $name = $1;
        if ($cols[$i] >= 0.1) {
          foreach my $gene (@genes) {
            if ($gene ne ''){
              $result{$gene}{$name}++;
            }                   #gene ne ''
          }                     #foreach
        }                       #found in this sample
      }                         #it is maf
    }                           #iterator
  }                             #else
}
close IN;
for (my $i = 0; $i <= $#name; $i++){
  if ($name[$i] =~ /(AC\d+)maf/){
     push(@samples, $1);
   }
}


print "gene";
if ($cohort == 1){
 foreach my $sample (@rectum) {
   print "\t$sample";
 }
 foreach my $sample (@ileum) {
   print "\t$sample";
 }
 print "\trectum\tileum\n";
} else {
  foreach my $sample (@samples) {
    print "\t$sample";
  }
  print "\n";
}

foreach my $gene (keys %result) {
   print "$gene";
   my $il = 0;
   my $rec = 0;
   foreach my $sample (keys %{$result{$gene}}) {
     if (exists($ileum{$sample})) {
       $il++;
     } elsif (exists($rectum{$sample})) {
       $rec++ unless ($sample eq 'AC442'); #do not count AC442
     }
   }
   if ($cohort == 1) {
     foreach my $sample (@rectum) {
       if ($result{$gene}{$sample} ne '') {
          print "\t$result{$gene}{$sample}";
       } else {
          print "\t0";
       }
     }
     foreach my $sample (@ileum) {
       if ($result{$gene}{$sample} ne ''){
          print "\t$result{$gene}{$sample}";
        } else {
          print "\t0";
        }
     }
     print "\t$rec\t$il\n";
   } else {
     foreach my $sample (@samples) {
       if ($result{$gene}{$sample} ne ''){
          print "\t$result{$gene}{$sample}";
       } else {
          print "\t0";
       }
     }
     print "\n";
   }
}

sub grepGene {
  my $line = shift;
  my @dummy;
  if ($line =~ /function=exonic\;functionalClass=([^\;]+)\;geneName=([\w\.\-\,\;\/]+?)\;[\w\.]+\=/) {
    my $g1 = $2;
    my $function = $1;
    if ($type =~ /exonic/){
      &splitGene($g1);
    } elsif ($type =~ /coding/) {
      if ($function =~ /^synonymous/ or $function =~ /^nonframeshift/) {
         return(@dummy);
      } else {
         &splitGene($g1);
      }
    }
  } elsif ($line =~ /function\=ncRNA\_exonic\;geneName=([\w\.\-\,\;\/]+?)\;[\w\.]+\=/) {
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /exonic/);
  } elsif ($line =~ /function\=ncRNA\_exonic\;geneName=([\w\.\-\,\;\/]+?)\t/) {
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /exonic/);
  } elsif ($line =~ /function=exonic\;splicing\;functionalClass=([^\;]+)\;geneName=([\w\.\-\,\;\(\)\_\:\+\>\/]+?)\;[\w\.]+\=/) {
    my $g1 = $2;
    my $function = $1;
    if ($type =~ /exonic/ or $type =~ /splicing/){
       &splitGene($g1);
    } elsif ($type =~ /coding/) {
      if ($function =~ /^synonymous/ or $function =~ /^nonframeshift/) {
        return(@dummy);
      } else {
        &splitGene($g1);
      }
    }
  } elsif ($line =~ /function\=(ncRNA\_)?intronic\;geneName=([\w\.\-\,\;\/]+?)\;[\w\.]+\=/) {
    my $g1 = $2;
    &splitGene($g1) if ($type =~ /intronic/);
  } elsif ($line =~ /function\=(ncRNA\_)?intronic\;geneName=([\w\.\-\,\;\/]+?)\t/) {
    my $g1 = $2;
    &splitGene($g1) if ($type =~ /intronic/);
  } elsif ($line =~ /function\=(ncRNA\_)?splicing\;geneName=([\w\.\-\,\;\(\)\_\:\+\>\/]+?)\;[\w\.]+\=/) {
    my $g1 = $2;
    &splitGene($g1) if ($type =~ /splicing/);
  } elsif ($line =~ /function\=(ncRNA\_)?splicing\;geneName=([\w\.\-\,\;\(\)\_\:\+\>\/]+?)\t/) {
    my $g1 = $2;
    &splitGene($g1) if ($type =~ /splicing/);
  } elsif ($line =~ /function\=(ncRNA\_)?UTR\d\;geneName=([\w\.\-\,\;\/]+?)\;[\w\.]+\=/) {
    my $g1 = $2;
    &splitGene($g1) if ($type =~ /UTR/);
  } elsif ($line =~ /function\=(ncRNA\_)?UTR\d\;geneName=([\w\.\-\,\;\/]+?)\t/){
    my $g1 = $2;
    &splitGene($g1) if ($type =~ /UTR/);
  } elsif ($line =~ /function\=UTR5\;UTR3\;geneName=([\w\.\-\,\;\/]+?)\;[\w\.]+\=/) {
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /UTR/);
  } elsif ($line =~ /function\=UTR5\;UTR3\;geneName=([\w\.\-\,\;\/]+?)\t/) {
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /UTR/);
  } elsif ($line =~ /function\=upstream\;geneName=([\w\.\-\,\;\/]+?)\;[\w\.]+\=/){
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /upstream/);
  } elsif ($line =~ /function\=upstream\;geneName=([\w\.\-\,\;\/]+?)\t/){
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /upstream/);
  } elsif ($line =~ /function\=downstream\;geneName=([\w\.\-\,\;\/]+?)\;[\w\.]+\=/){
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /downstream/);
  } elsif ($line =~ /function\=downstream\;geneName=([\w\.\-\,\;\/]+?)\t/){
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /downstream/);
  } elsif ($line =~ /function\=upstream\;downstream\;geneName=([\w\.\-\,\;\/]+?)\;[\w\.]+\=/){
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /upstream/ or $type =~ /downstream/);
  } elsif ($line =~ /function\=upstream\;downstream\;geneName=([\w\.\-\,\;\/]+?)\t/) {
    my $g1 = $1;
    &splitGene($g1) if ($type =~ /upstream/ or $type =~ /downstream/);
  } elsif ($line =~ /function\=intergenic\;geneName\=(.+?)\(dist\=\w+\)\,(.+?)\(dist\=\w+\)/) {
     my $g1 = $1;
     my $g2 = $2;
     &splitGene($g1) if ($type =~ /intergenic/ and $g1 ne 'NONE');
     &splitGene($g2) if ($type =~ /intergenic/ and $g2 ne 'NONE');
  } else {
    print STDERR "error:\t\t$_\n";
  }
}


sub splitGene {
  my $genes = shift;
  $genes =~ s/\(.+?\)//g;    #remove the splicing parentheses
  my @genes = split (/[\;\,]/, $genes);
  my %tmp;
  foreach my $gene (@genes) {
    $tmp{$gene} = "";
  }
  @genes = keys %tmp;
  return(@genes);
}
