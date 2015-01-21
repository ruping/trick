use strict;
use Getopt::Long;
use List::Util qw(min max);

my $mutation;
my $type = 'snv';
my $normal;
my $task; #arj cohort or other?
my $prefix;
my $tmpdir = './';
my $clinical = 0;

GetOptions (
           "mutation|m=s"   => \$mutation,          #filename of the mutation table
           "type|t=s"       => \$type,             #snv or indel
           "normal|n=s"     => \$normal,           #comma seperated id of normal samples
           "task|k=s"       => \$task,             #task type
           "prefix|p=s"     => \$prefix,
           "tmpdir|y=s"     => \$tmpdir,
           "clinical=i"     => \$clinical,
           "help|h"         => sub{
                               print "usage: $0 produce andrea wanted mutation table\n\nOptions:\n\t--mutation\tthe filename of mutation table\n";
                               print "\t--type\t\tthe type of variants, snv or indel\n";
                               print "\t--normal\tcomma seperated id of normal samples\n";
                               print "\t--prefix\tthe prefix of samples' names\n";
                               print "\t--task\t\tthe task, such as tcga or rnaediting\n";
                               print "\t--tmpdir\tthe temporary dir to write tmp files\n";
                               print "\t--clinical\twhether (1) or not (0 default) only print clinical vars\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );


my @type1 = qw(AC2maf AC3maf AC4maf AC53maf AC54maf AC55maf AC56maf AC57maf AC565maf AC566maf AC567maf) if ($task eq 'arj');
my %type1;
foreach my $sample (@type1){
  $type1{$sample} = 0;
}
my @type2 = qw(AC582maf AC583maf) if ($task eq 'arj');
my %type2;
foreach my $sample (@type2){
  $type2{$sample} = 0;
}

open IN, "$mutation";
my @colnames;
while ( <IN> ) {
  chomp;
  if ($_ =~ /^[\#]?[cC]hr\t/) {
    @colnames = split /\t/;
    for(my $c = 0; $c <= $#colnames; $c++){
      if ($c == 0) {
        print "$colnames[$c]";
      } elsif ($colnames[$c] eq 'id'){
        print "\tlink\tid";
      } elsif ($colnames[$c] =~ /function/){
        print "\tgeneName\tgeneLoc\tfunctionalClass\tAAChange";
      } elsif ($colnames[$c] eq 'clinical'){
        print "\tpopFreq\t\tClinChanel\tClinAllele\tClinVariantDisease";
      } else {
        print "\t$colnames[$c]";
      }
    }
    print "\tcloneType\n" if ($task eq 'arj');
    print "\n" if ($task eq 'cohort' or $task eq 'other');
  } else {
    my @cols = split /\t/;
    my @printcols;
    my $print;
    my $onoff = 1;
    my $chr;
    my $pos;
    my $ucsc;
    my $nt1 = 0;
    my $nt2 = 0;
    my $nt1nz = 0;
    my $nt2nz = 0;
    my $primary = 0;
    my $cloneType = "somatic";
    for (my $i = 0; $i <= $#colnames; $i++) {
      if ($colnames[$i] eq 'chr') {
        $chr = $cols[$i];
        $print = $chr;
        push (@printcols, $print);

      } elsif ($colnames[$i] eq 'pos') {
        $pos = $cols[$i];
        $ucsc = "\=HYPERLINK\(\"http\:\/\/genome\.ucsc\.edu\/cgi\-bin/hgTracks\?db\=hg19\&position\=chr$chr\%3A$pos\-$pos\", \"UCSC\"\)";
        $print = $pos."\t".$ucsc;
        push (@printcols, $print);

      } elsif ($colnames[$i] eq 'id' or $colnames[$i] eq 'ref' or $colnames[$i] eq 'alt') {
        $print = $cols[$i];
        push (@printcols, $print);

      } elsif ($colnames[$i] =~ /function/) { #now it is function, gene names need to be extracted
        my $functions = &splitFunction($cols[$i]);
        $print = "$functions->{'geneName'}\t$functions->{'loc'}\t$functions->{'functionClass'}\t$functions->{'AAChange'}";
        push (@printcols, $print);

      } elsif ($colnames[$i] eq 'clinical') {
        my $clins = &splitClinical($cols[$i]);
        $print = "$clins->{'freq'}\t$clins->{'clinVarChanel'}\t$clins->{'clinAllele'}\t$clins->{'clinVarDisease'}";
        push (@printcols, $print);
        if ($printcols[2] eq '.' and $clins->{'RS'} ne 'NA'){
          $printcols[2] = $clins->{'RS'};
        }
        if ($clinical == 1 and $clins->{'clinVarChanel'} eq 'NA'){
          $onoff = 0;
        }

      } elsif ( $colnames[$i] =~ /$prefix\d+maf/ ) {
        print "\t$cols[$i]";
        if (exists ($type1{$colnames[$i]})) {
          $nt1++ if ($cols[$i] >= 0.1);
          $nt1nz++ if ($cols[$i] > 0);
        } elsif (exists ($type2{$colnames[$i]})) {
          $nt2++ if ($cols[$i] >= 0.1);
          $nt2nz++ if ($cols[$i] > 0);
        }
        if ( $colnames[$i] eq 'AC565maf' and $cols[$i] == 0 ) {
          $primary += 0.5;
        }

      }

      #elsif ( $colnames[$i] =~ /^AC\d+d$/ ) {
      #  if ( $colnames[$i] eq 'AC565d' and $cols[$i] >= 5 ) {
      #    $primary += 0.5;
      #  }
      #  $print = $cols[$i];
      #  push (@printcols, $print);
      #
      #}

      else {
        $print = $cols[$i];
        push (@printcols, $print);
      }
    }
    #if (($nt1 >= 1 and $nt1nz >= 5) and ($nt2 > 0 and $nt2nz > 1)) {
    #  $cloneType .= '.both'
    #} elsif (($nt1 <= 1 and $nt1nz < 3) and ($nt2 > 0 and $nt2nz > 1)){
    #  $cloneType .= '.type2';
    #} elsif (($nt2 == 0 and $nt2nz == 0) and $nt1 >= 2){
    #  $cloneType .= '.type1';
    #}
    #if ($primary == 1){
    #  $cloneType .= '.metastasis';
    #}
    #if ($task eq 'arj'){
    #  $print = $cloneType;
    #  push (@printcols, $print);
    #}

    if ($onoff == 1){
      printf("%s\n", join("\t",@printcols));
    }
  }
}
close IN;



sub splitFunction {
  my $func = shift;
  my %func = (
              'loc' => 'NA',
              'functionClass' => 'NA',
              'geneName' => 'NA',
              'AAchange' => 'NA',
              );
  if ($func =~ /function=([^\;]+)/){
    $func{'loc'} = $1;
  }
  if ($func =~ /functionalClass=([^\;]+)/){
    $func{'functionClass'} = $1;
  }
  if ($func =~ /geneName=([^\;]+)/){
    $func{'geneName'} = $1;
  }
  if ($func =~ /AAChange=([^\;]+)/){
    $func{'AAChange'} = $1;
  }
  return(\%func);
}


sub splitClinical {
  my $clin = shift;
  my %clin = (
              'RS' => 'NA',
              'freq' => 'NA',
              'clinAllele' => 'NA',
              'clinVarChanel' => 'NA',
              'clinVarDisease' => 'NA',
             );
  if ($clin =~ /RS=([^\;]+)/){
    $clin{'RS'} = 'rs'.$1;
  }
  if ($clin =~ /CLNALLE=([01])/){
    $clin{'clinAllele'} = $1;
  }
  if ($clin =~ /CLNSRC=([^\;]+)/){
    $clin{'clinVarChanel'} = $1;
  }
  if ($clin =~ /CLNDBN=([^\;]+)/){
    $clin{'clinVarDisease'} = $1;
  }
  if ($clin =~ /CAF\=\[([\d\.\,]+)\]/){
    my @freqs = split (/\,/, $1);
    shift @freqs;
    $clin{'freq'} = max(@freqs);
  }
  return(\%clin);
}
