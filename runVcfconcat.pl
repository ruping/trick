use strict;
use File::Glob ':glob';
use File::Basename;
use Data::Dumper;

my $dir = shift;
my $outdir = shift;
my $type = shift;
my @vcfs;
if ($type eq 'muTect') {
  @vcfs = bsd_glob("$dir/*.vcf");
} elsif ($type eq 'strelka') {
  @vcfs = bsd_glob("$dir/*/results/all.somatic.indels.vcf");
}
#print STDERR Dumper(\@vcfs);


my %vcfs;
my $outvcfname = '';
foreach my $vcf (@vcfs) {
  my $dirname = dirname($vcf);
  my $basename = basename($vcf);
  my $sample;
  my $chr;
  my $fname;
  if ($type eq 'muTect' and $basename =~ /^([A-Za-z0-9\-\_]+)\.([A-Za-z0-9]+)\.(.+?)$/) {
    $sample = $1;
    $chr = $2;
    $fname = $3;
    if ($outvcfname eq '') {
      $outvcfname = $sample.'.'.$fname;
    }
    $vcfs{$chr} = $vcf;
  } elsif ($type eq 'strelka' and $dirname =~ /\/([A-Za-z0-9\-\_]+)(\/)+([A-Za-z0-9]+)(\/)+results\//) {
    $sample = $1;
    $chr = $3;
    $fname = 'genome.somatic.indel.vcf';
    if ($outvcfname eq '') {
      $outvcfname = $sample.'.'.$fname;
    }
    $vcfs{$chr} = $vcf;
  } else {
    print "strang $vcf\n";
    exit 22;
  }
}

print STDERR Dumper(\%vcfs);


=pod
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
