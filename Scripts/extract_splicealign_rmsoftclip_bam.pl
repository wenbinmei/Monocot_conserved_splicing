#!/usr/bin/perl -w 
use strict;
use warnings;
use Bio::DB::Sam;

### This scripts extract the splice alignment from the bam file

my $usage = "$0 <File>";
die $usage unless $#ARGV == 1;
my ($File, $output) = @ARGV;

my $out = Bio::DB::Bam->open($output,'w')
	or die "Could not open BAM file for writing: $!";

### read the bam alignment
my $bam = Bio::DB::Bam->open("$File");
my $header       = $bam->header;
my $target_names = $header->target_name;
$out->header_write($header);

while (my $align = $bam->read1) {
        my $chr       = $target_names->[$align->tid];
        my $start     = $align->pos+1;
        my $end       = $align->calend;
        my $cigar     = $align->cigar_str;
      
        if ( ($cigar =~ /N/) && ($cigar !~ /S/) ) { ### extract the reads with splice alignment
             $out->write1($align);
        }
}
  
