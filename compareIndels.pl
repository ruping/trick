use strict;
use Data::Dumper;

my $tumorSam = shift;
my $normalSam = shift;

open TS, "$tumorSam";
my %tumorIndels;         #save the tumor indel data
while ( <TS> ){
  next if /^@/;
  chomp;
  my @cols = split /\t/;
  $cols[0] =~ /id\_(\d+)$/;
  my $id = $1;

  my $chr = $cols[2];
  $chr = 'chr'.$chr if ($chr !~ /^[cC]hr/);
  my $cigar = $cols[5];
  my @cigarData;
  my $alignmentStart = $cols[3];

  while ($cigar =~ /(\d+)([MIDNSHPX=])/g){
    my $clen = $1;
    my $ctype = $2;
    push(@cigarData, {"Length"=>$clen, "Type"=>$ctype});
  }

  my $cigarEnd;
  my @blockLengths;
  my @blockStarts;
  my %insertions;       # for insertions
  my %deletions;        # for deletions
  my $softClip = 0;     # for soft clipping
  push(@blockStarts, 0);
  &ParseCigar(\@cigarData, \@blockStarts, \@blockLengths, \$cigarEnd, \%insertions, \%deletions, \$softClip);
  my $alignmentEnd = $alignmentStart + $cigarEnd - 1;

  #debug###########################
  #print STDERR Dumper(\@cigarData);
  #print STDERR Dumper(\@blockLengths);
  #print STDERR Dumper(\@blockStarts);
  #print STDERR Dumper(\%insertions);
  #print STDERR Dumper(\%deletions);
  #print STDERR "$cigarEnd\t$alignmentEnd\t$softClip\n";

  my $indelType;
  my $indelSite;
  my $indelLength = 0;
  foreach my $insStart (sort {$a <=> $b} keys %insertions) {
    my $insPos = $alignmentStart + $insStart - 1;
    my $insLen = $insertions{$insStart};
    if ($insLen > $indelLength){
      $indelType = 'I';
      $indelSite = $insPos;
      $indelLength = $insLen;
    }
  }

  foreach my $delStart (sort {$a <=> $b} keys %deletions) {
    my $delPos = $alignmentStart + $delStart - 1;
    my $delLen = $deletions{$delStart};
    if ($delLen > $indelLength){
      $indelType = 'D';
      $indelSite = $delPos;
      $indelLength = $delLen;
    }
  }

  #print "$id\t$chr\t$alignmentStart\t$alignmentEnd\t$cigar\t$indelType\t$indelSite\t$indelLength\n";
  if ($indelType ne '') {
    push (@{$tumorIndels{$id}}, {'chr'=>$chr, 'site'=>$indelSite, 'type'=>$indelType, 'length'=>$indelLength, 'germline'=>0});
  }

}
close TS;


#look at the normal sample#####################################################################################

open NS, "$normalSam";
while ( <NS> ){
  next if /^@/;
  chomp;
  my @cols = split /\t/;
  $cols[0] =~ /id\_(\d+)$/;
  my $id = $1;

  my $chr = $cols[2];
  $chr = 'chr'.$chr if ($chr !~ /^[cC]hr/);
  my $cigar = $cols[5];
  my @cigarData;
  my $alignmentStart = $cols[3];

  while ($cigar =~ /(\d+)([MIDNSHPX=])/g){
    my $clen = $1;
    my $ctype = $2;
    push(@cigarData, {"Length"=>$clen, "Type"=>$ctype});
  }

  my $cigarEnd;
  my @blockLengths;
  my @blockStarts;
  my %insertions;       # for insertions
  my %deletions;        # for deletions
  my $softClip = 0;     # for soft clipping
  push(@blockStarts, 0);
  &ParseCigar(\@cigarData, \@blockStarts, \@blockLengths, \$cigarEnd, \%insertions, \%deletions, \$softClip);
  my $alignmentEnd = $alignmentStart + $cigarEnd - 1;

  my $indelType;
  my $indelSite;
  my $indelLength = 0;
  foreach my $insStart (sort {$a <=> $b} keys %insertions) {
    my $insPos = $alignmentStart + $insStart - 1;
    my $insLen = $insertions{$insStart};
    if ($insLen > $indelLength){
      $indelType = 'I';
      $indelSite = $insPos;
      $indelLength = $insLen;
    }
  }

  foreach my $delStart (sort {$a <=> $b} keys %deletions) {
    my $delPos = $alignmentStart + $delStart - 1;
    my $delLen = $deletions{$delStart};
    if ($delLen > $indelLength){
      $indelType = 'D';
      $indelSite = $delPos;
      $indelLength = $delLen;
    }
  }

  if (exists ($tumorIndels{$id})){
    foreach my $tumorIndel (@{$tumorIndels{$id}}) {
      my $tichr = $tumorIndel->{'chr'};
      my $titype = $tumorIndel->{'type'};
      my $tisite = $tumorIndel->{'site'};
      my $tilength = $tumorIndel->{'length'};
      if ( ($chr eq $tichr) and ($indelType eq $titype ) and (abs($indelSite-$tisite) <= 10) and ($indelLength eq $tilength) ) {
        $tumorIndel->{'germline'} = 1;
      }
    }
  }
}
close NS;
#look at the normal sample#####################################################################################

foreach my $id (sort {$a <=> $b} keys %tumorIndels) {
  foreach my $tumorIndel (@{$tumorIndels{$id}}) {
    my $tichr = $tumorIndel->{'chr'};
    my $titype = $tumorIndel->{'type'};
    my $tisite = $tumorIndel->{'site'};
    my $tilength = $tumorIndel->{'length'};
    my $germline = $tumorIndel->{'germline'};
    if ($germline == 0){
      printf("%s\n",join("\t", $id, $tichr, $titype, $tisite, $tilength, $germline));
    }
  }
}

sub ParseCigar {  #process cigar string

  my ($cigar, $blockStarts, $blockLengths, $cigarEnd, $insertions, $deletions, $softClip) = @_;

  my $currPosition = 0;
  my $blockLength  = 0;

  my $insertSize;
  my $insertPos;
  my $deleteSize;
  my $deletePos;

  #  Rip through the CIGAR ops and figure out if there is more than one block for this alignment
  for (my $cigItr = 0; $cigItr <= $#{$cigar}; $cigItr++) {
    if ($$cigar[$cigItr]->{"Type"} eq 'M') { # matching
      $blockLength += $$cigar[$cigItr]->{"Length"};
      $currPosition += $$cigar[$cigItr]->{"Length"};
    } elsif ($$cigar[$cigItr]->{"Type"} eq 'I') { #  insertion
      $insertSize = $$cigar[$cigItr]->{"Length"};
      $insertPos = $currPosition + 1;
      $insertions->{$insertPos} = $insertSize;
    } elsif ($$cigar[$cigItr]->{"Type"} eq 'S') { # soft-clipping
      if ($currPosition == 0) { #only take action for the beginning clipping
        $$softClip = $$cigar[$cigItr]->{"Length"};
      }
    } elsif ($$cigar[$cigItr]->{"Type"} eq 'D') { # deletion
      $deleteSize = $$cigar[$cigItr]->{"Length"};
      $deletePos = $currPosition + 1;
      $deletions->{$deletePos} = $deleteSize;
      $blockLength  += $$cigar[$cigItr]->{"Length"};
      $currPosition += $$cigar[$cigItr]->{"Length"};
    } elsif ($$cigar[$cigItr]->{"Type"} eq 'P') {
      #do nothing
    } elsif ($$cigar[$cigItr]->{"Type"} eq 'N') { # skipped region
      push(@{$blockStarts}, $currPosition + $$cigar[$cigItr]->{"Length"});
      push(@{$blockLengths}, $blockLength);
      $currPosition += $$cigar[$cigItr]->{"Length"};
      $blockLength = 0;         # a new block
    } elsif ($$cigar[$cigItr]->{"Type"} eq 'H') {
      #do nothing
    } else {
      print STDERR "ERROR: Invalid Cigar op type\n"; # shouldn't get here
      exit 22;
    }
  }
  # add the kast block and set the
  # alignment end (i.e., relative to the start)
  push(@{$blockLengths}, $blockLength);
  $$cigarEnd = $currPosition;

}
