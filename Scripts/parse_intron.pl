#!/usr/bin/perl -w
use strict;
use warnings;

### >evm_15.TU.AmTr_v1.0_scaffold00033.1.evm_15.model.AmTr_v1.0_scaffold00033.1.intron1/19 AmTr_v1.0_scaffold00033:151567..152658
### >PASA_cluster_1.align_id:1721244|asmbl_1.intron1/8 chr1:9193..7903

my $usage = "$0 <File>";

die $usage unless $#ARGV == 0;

my ($File) = @ARGV;

print "Intron_Name\tChr\tStart\tEnd\tSize\n";

open (FILE, $File) or die "Can not open input file: $File\n";
while (<FILE>) {
      chomp;
      my @line = split (/\s+/, $_);
      my $Intron_Name = $line[0];
      $Intron_Name =~ s/>//g;  ### remove >
      my @temp = split (/\:/, $line[1]);
      my $chr = $temp[0];
      my $location = $temp[1];
      $location =~ s/\../_/g;  ### convert .. to _ which is easy to split
      #print "$temp[1]\n"; 
      my @array = split (/\_/, $location);
      #print "$array[0]\n";
      if ($array[0] > $array[1]) { #### this intron is on the negative strand 
          my $intron_start = $array[1] + 1;
          my $intron_end = $array[0] - 1;
          my $distance = abs($intron_end - $intron_start) + 1;
          print "$Intron_Name\t$chr\t$intron_start\t$intron_end\t$distance\n";
      }
      else {  ### intron on the positive strand
          my $intron_start = $array[0] + 1;
          my $intron_end = $array[1] - 1;    ### coordinates are 1-based, and specify the exon coordinates surrounding the intron, with the first coordinate being from the donor exon and the second one being from the acceptor exon.
          my $distance = abs($intron_end - $intron_start) + 1;
          print "$Intron_Name\t$chr\t$intron_start\t$intron_end\t$distance\n";
      }
}
close FILE;
