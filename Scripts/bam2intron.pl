#!/usr/bin/perl -w 
###############################################
### bam2intron.pl
### Wenbin Mei 
### Barbazuk lab
### The script take the bam file and output the junction file with counts, the start and the end of the coordinate is the first base and last base of the intron.
###############################################


use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Sam;

my ($infile, $source, $output, $help);
GetOptions("input|i=s" => \$infile, "source|s=s" => \$source, "output|o=s" => \$output, "help|h" => \$help);

if ($help) {
	print <<EOF;
Usae: perl bam2intron.pl -i [BAM file] -s [source string]
		Options
		--input|i: bam file
		--source|s: string to add in the header
                --output|o: output file name
		--help: print help information
EOF
exit;
}

#print header:
my (%intron, @line);

### read the bam alignment
my $bam = Bio::DB::Bam->open("$infile");
my $header       = $bam->header;
my $target_names = $header->target_name;
while (my $align = $bam->read1) {
	my $chr	      = $target_names->[$align->tid];
        my $start     = $align->pos+1;
        my $end       = $align->calend;
        my $cigar     = $align->cigar_str;

        if ($cigar =~ /N/) { ### extract the reads with splice alignment
                  my @large_gap = split(/N/,$cigar);
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

open (OUTPUT, ">", "$output"); 

#print header:
print OUTPUT "chr\tintronStart\tintronEnd\t$source\n";

#output for each chr:
foreach my $eachintron (sort {$a cmp $b} keys %intron) {
         print OUTPUT "$eachintron\t$intron{$eachintron}\n";
}

system("sort -k1,1 -k2,2n -k3,3n $output > $output.'sorted'");
system("mv $output.'sorted' $output.'junction' ");
unlink $output.'sorted';

