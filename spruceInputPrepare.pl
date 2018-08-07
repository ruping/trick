use strict;
use List::Util qw[max];

open IN, shift;
my $header;
my $m;
my $n;
my @allsampindex;
my @allmutindex;
my %segs;
my %data;
while ( <IN> ){
  chomp;
  if (/^sample\_index/){   #header
    $header = $_;
    $header =~ s/\tseg$//;
    $header = '#'.$header;
  } else {
    my ($sample_index, $sample_label, $character_index, $character_label, $vaf_lb, $vaf_mean, $vaf_ub, $x1, $y1, $mu1, $x2, $y2, $mu2, $seg) = split /\t/;
    my %tmp = (
               'sample_label' => $sample_label,
               'character_label' => $character_label,
               'vaf_lb' => $vaf_lb,
               'vaf_mean' => $vaf_mean,
               'vaf_ub' => $vaf_ub,
               'x1' => $x1,
               'y1' => $y1,
               'mu1' => $mu1,
               'x2' => $x2,
               'y2' => $y2,
               'mu2' => $mu2,
               'seg' => $seg
               );
    $data{$character_index}{$sample_index} = \%tmp;
    push(@allsampindex, $sample_index);
    push(@allmutindex, $character_index);

    if ($seg != 0) {
      $segs{$sample_index}{$seg} .= $character_index." ";
    }
  }
}
close IN;

$m = max(@allsampindex)+1;
$n = max(@allmutindex)+1;
print "$m #m\n$n #n\n";

my $maxnanbs = 2;   #max number of states for all mutations
my @printinfo;
foreach my $character_index (sort {$a <=> $b} keys %data) {  #each mutation

  my %nanbs;  #store all type of nanb for this mutation
  foreach my $sample_index (sort {$a <=> $b} keys %{$data{$character_index}}) {  #each sample
    my $info = $data{$character_index}{$sample_index};  #this is a hash table
    if ($info->{'x1'} ne 'NA'){
      $nanbs{$info->{'x1'}.$info->{'y1'}}++;
    }
    if ($info->{'x2'} ne 'NA'){
      $nanbs{$info->{'x2'}.$info->{'y2'}}++;
    }
  }
  if ( scalar(keys %nanbs) > $maxnanbs ){
    $maxnanbs = scalar(keys %nanbs);    #max number of states for all mutations
  }

  foreach my $sample_index (sort {$a <=> $b} keys %{$data{$character_index}}) {  #print data
    my $info = $data{$character_index}{$sample_index};  #this is a hash table
    my $printinfo = join("\t", $sample_index, $info->{'sample_label'}, $character_index, $info->{'character_label'},
                         $info->{'vaf_lb'}, $info->{'vaf_mean'}, $info->{'vaf_ub'});

    my $nanb1 = 'SRP';
    $nanb1 = $info->{'x1'}.$info->{'y1'} if ($info->{'x1'} ne 'NA');
    my $nanb2 = 'SRP';
    $nanb2 = $info->{'x2'}.$info->{'y2'} if ($info->{'x2'} ne 'NA');
    my $printadd = '';
    foreach my $nanb (sort keys %nanbs) {
      if ($nanb1 ne 'SRP' and $nanb1 eq $nanb) {
        $nanb =~ /^(\d)(\d)$/;
        $printadd .= "\t".join("\t",$1,$2,$info->{'mu1'});
      } elsif ($nanb2 ne 'SRP' and $nanb2 eq $nanb) {
        $nanb =~ /^(\d)(\d)$/;
        $printadd .= "\t".join("\t",$1,$2,$info->{'mu2'});
      } else {
        $nanb =~ /^(\d)(\d)$/;
        $printadd .= "\t".join("\t",$1,$2,0);
      }
    }
    $printinfo .= $printadd;
    push(@printinfo, $printinfo);
  }

}

#print header
printf("%s\n",join("", $header, "\tx\ty\tmu"x($maxnanbs - 2)));
#print content
foreach my $printinfo (@printinfo){
  print "$printinfo\n";
}


foreach my $sample (sort {$a <=> $b} keys %segs){
  my $linked;
  foreach my $seg (sort {$a <=> $b} keys %{$segs{$sample}}) {
    my $cl= $segs{$sample}{$seg};
    $cl =~ s/\s$//;
    $linked .= "$cl\t"
  }
  $linked =~ s/\t$//;
  print STDERR "$linked\n";
}

exit 0;
