#!/usr/bin/perl -w 
# sam2intron.pl
# Sanzhen Liu
# 7/18/2011

use strict;
use warnings;
use Getopt::Long;

my ($infile, $source, $help);
GetOptions("input|i=s" => \$infile, "source|s=s" => \$source, "help|h" => \$help);

if ($help) {
	print <<EOF;
Usae: perl sam2intron.pl -i [SAM file] -s [source string]
		Options
		--input|i: SAM file
		--source|s: string to add in the header
		--help: print help information
EOF
exit;
}

#print header:
print "chr\tintronStart\tintronEnd\t$source\n";
my (%intron, @line, $chr, $start);
open(IN,$infile) || die;
while (<IN>) {
	if (!/^@/) {
		chomp;
		@line = split(/\t/,$_);
		$chr = $line[2];
		$start = $line[3];
		my $align = $line[5]; # alignment summary
		if ($align =~ /N/) {
			my @large_gap = split(/N/,$align);
			my $pos = $start;
			for (my $i=0; $i<$#large_gap; $i++) {
				$large_gap[$i] =~ /(\d+)$/;
				my $intron_size = $1;
				my @typecount = split(/[MIDSHP]/,$large_gap[$i]); # array of matched, gap bases and so on
				$large_gap[$i] =~ s/[0-9]//g; # remove those number
				my @type = split(//,$large_gap[$i]);
				for (my $i=0; $i<=$#type; $i++) {
					# gap, deletion:
					if ($type[$i] =~ /[MD]/) {  # deletion
						$pos+=$typecount[$i];
					}
				} # end FOR
				my $intron_start = $pos;
				my $intron_end = $pos+$intron_size-1;
				$pos = $intron_end + 1;
				my $intron_id = $chr."\t".$intron_start."\t".$intron_end;
				$intron{$intron_id}++;
			}
		}
	}
} # end of while <IN>
close IN;

# output for each chr:
foreach my $eachintron (sort {$a cmp $b} keys %intron) {
	print "$eachintron\t$intron{$eachintron}\n";
}


#### subroutines ####
sub compare {
	# input two numbers: 1. to-be checked number 2. bit number
	my $out=0;
	my ($in,$bit) = @_;
	if ($in>=$bit) {
		$out=1;
	}
}

sub decode {
	# two factors are input:
	# first one is the number to-be judged
	# second one is the code to-be known, e.g. duplication or seq orientation
	my ($num,$item) = @_;
	my $value=0;

	my $yesno = compare($num,1024);
	if ($yesno) {$num-=1024; $yesno=0; if ($item==1024) { $value=1; }} # 0x400: PCR or optical duplicate
	
	$yesno = compare($num,512);
	if ($yesno) {$num-=512; $yesno=0; if ($item==512) { $value=1; }} # 0x200: not passing quality controls
	
	$yesno = compare($num,256);
	if ($yesno) {$num-=256; $yesno=0; if ($item==256) { $value=1; }} # 0x100: secondary alignment
	
	$yesno = compare($num,128);
	if ($yesno) {$num-=128; $yesno=0; if ($item==128) { $value=1; }} # 0x80: the last fragment in the template
	
	$yesno = compare($num,64);
	if ($yesno) {$num-=64; $yesno=0; if ($item==64) { $value=1; }} # 0x40: the first fragment in the template
	
	$yesno = compare($num,32);
	if ($yesno) {$num-=32; $yesno=0; if ($item==32) { $value=1; }} # 0x20: SEQ of the next fragment in the template being reversed

	$yesno = compare($num,16);
	if ($yesno) {$num-=16; $yesno=0; if ($item==16) { $value=1; }} # 0x10: SEQ being reversed complemented

	$yesno = compare($num,8);
	if ($yesno) {$num-=8; $yesno=0; if ($item==8) { $value=1; }} # 0x8: next fragment in the template unmapped

	$yesno = compare($num,4);
	if ($yesno) {$num-=4; $yesno=0; if ($item==4) { $value=1; }} # 0x4: fragment unmapped

	$yesno = compare($num,2);
	if ($yesno) {$num-=2; $yesno=0; if ($item==2) { $value=1; }} # 0x2: each fragment properly aligned according to the aligner
	
	$yesno = compare($num,1);
	if ($yesno) {$yesno=0; if ($item==1) { $value=1; }} # 0x1: template having multiple fragments in sequencing
	
	return $value;
}

# returns the minimum value from an array of numbers
sub min {
	my @sorted = sort {$a <=> $b} @_;
	return shift(@sorted);
} # End of sub min

# returns the maximum value from an array of numbers
sub max {
	my @sorted = sort {$a <=> $b} @_;
	return pop(@sorted);
} # End of sub max



