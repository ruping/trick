use strict;
use File::Glob ':glob';
use File::Basename;
use Data::Dumper;

my $dir = shift;
my $outdir = shift;
my $type = shift;
my @vcfs;
if ($type eq 'muTect'){
  @vcfs = bsd_glob("$dir/*.vcf");
} elsif ($type eq ''){
  @vcfs = bsd_glob("$dir/*/results/all.somatic.indels.vcf");
}
print STDERR Dumper(\@vcfs);

=pod
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
=cut
