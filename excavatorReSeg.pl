use strict;
use Data::Dumper;
use File::Glob ':glob';


my $dir = shift;
my $prefix = shift;

my @files = bsd_glob("$dir/$prefix*/FastCallResults_*.txt");
#print STDERR Dumper (\@files);

print "#FID\tCHR\tBP1\tBP2\tTYPE\n";

foreach my $file (@files) {
  my $sample;
  if ($file =~ /$dir\/($prefix\d+)\/FastCallResults_$prefix\d+\.txt/){
    $sample = $1;
  } else {
    print STDERR "$file is strange!!!\n";
    next;
  }
  #print STDERR "$sample\n";

  open IN, "$file";
  while ( <IN> ){
    chomp;
    next if /^Chromosome/;
    my ($chr, $start, $end, $logR, $CN, $Call, $ProbCall) = split /\t/;
    next if ($chr =~ /[XY]/);
    printf("%s\n", join("\t",$sample,$chr,$start,$end,$call));
  }
  close IN;
}

exit 0;
