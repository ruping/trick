use strict;
use File::Glob ':glob';
use File::Basename;

my $dir = shift;
my $outdir = shift;
my @vcfs = bsd_glob("$dir/*.vcf");

my %vcfs;
my $outvcfname = '';
foreach my $vcf (@vcfs) {
  my $basename = basename($vcf);
  if ($basename =~ /^([A-Za-z0-9\-\_]+)\.([A-Za-z0-9]+)\.(.+?)$/) {
    my $sample = $1;
    my $chr = $2;
    my $fname = $3;
    if ($outvcfname eq '') {
      $outvcfname = $sample.'.'.$fname;
    }
    $vcfs{$chr} = $vcf;
  } else {
    print "strang $vcf\n";
    exit 22;
  }
}

my $vcfs = '';
foreach my $chr (sort {$a cmp $b} keys %vcfs) {
  $vcfs .= "$vcfs{$chr} ";
}

my $cmd = "vcf-concat $vcfs >$outdir/$outvcfname";
unless (-s "$outdir/$outvcfname"){
  RunCommand($cmd, 0, 0);
}

sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}
