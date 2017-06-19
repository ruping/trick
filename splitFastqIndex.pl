use strict;
use File::Basename;
use Data::Dumper;

my $fastq = shift;
my $indexes = shift;
my $samples = shift;

my @indexes = split(',',$indexes);
my @samples = split(',',$samples);
my %indexes;
for (my $i = 0; $i <= $#indexes; $i++) {
  $indexes{$indexes[$i]} = $samples[$i];
}
#print STDERR Dumper(\%indexes);

my $bname = basename($fastq);
$bname =~ /Undetermined(_.+?\.fastq)\.gz$/;
my $suffix = $1;

my %fhs; #output
open IN, "gzip -dc $fastq |";
while ( <IN> ) {

  chomp;
  if ($_ =~ /^\@\S+\s+\S+\:([A-Z]+)$/) {
    my $cindex = $1;
    my $csample = $indexes{$cindex};
    if ($csample eq ''){$csample = "unknown"}
    my $fh = $csample;
    unless (-e "./$csample"."$suffix") {
      open ( my $fh, ">>", "./$csample"."$suffix" )  || die $!;
      $fhs{$csample} = $fh;
    }
    print {$fhs{$csample}} "$_\n";
    my $seq = <IN>;
    print {$fhs{$csample}} "$seq";
    my $third = <IN>;
    print {$fhs{$csample}} "$third";
    my $qual = <IN>;
    print {$fhs{$csample}} "$qual";
  }

}
close IN;

exit 0;
