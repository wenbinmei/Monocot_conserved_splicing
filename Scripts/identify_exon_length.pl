### this script take two files PASA_AS_File and Conserved_Event_File plus the species name
### and identify the span exon length (not the length of exon skipped) for skipped exon events which exist in the conserved AS or NOT

#Final_PASA_OUTPUT_pasa_wmei_B73_March18th_2015.indiv_splice_labels_and_coords.dat
#chr1    skipped_exon    -       GRMZM2G306216   asmbl_54        147671  146634  asmbl_51        147543  147423

#B73_event_ExonS_conserved_fiveASevents_10species_20Jan2016_table_RemoveUniqueConservedB73W22.txt
#GRMZM2G061988:ExonS_chr1_-_140050454_140048199_NA_NA
#GRMZM2G045314:ExonS_chr1_-_258005410_258002531_NA_NA

#!/usr/bin/perl -w
use strict;
use warnings;

my $usage = "$0 <PASA_AS_File> <Conserved_Event_File> <species name>";
die $usage unless $#ARGV == 2;
my ($AS_file, $conserved_AS, $species_name) = @ARGV;
my %conserved_AS;
my %output;

open (CONSERVED, $conserved_AS) or die "Can not open input file: $conserved_AS\n";
while (<CONSERVED>) {
		chomp;
		my ($gene, $info) = split (/\:/, $_);
		my @line = split (/\_/, $info);
		my $start = $line[-4];
		my $end = $line[-3];
		#print "$chr\t$strand\t$start\t$end\n";
		if (not defined $conserved_AS{$gene}{$start}{$end}) {
			$conserved_AS{$gene}{$start}{$end} = $gene."_".$start."_".$end;
		}
}
close CONSERVED;

print "Exon_span_length\tSpecies\tConserved\n";

open (AS, $AS_file) or die "Can not open input file: $AS_file\n";
while (<AS>) {
		next if ($_ !~ /skipped_exon/);
		my ($Chr, $AS_type, $strand, $Gene_ID, $Isoform_A, $A_start, $A_end, $Isoform_B, $B_start, $B_end) = split (/\s+/, $_);
                my $exon_span_size = abs($A_start - $A_end) + 1; ### because the coordinate is based on exon coordinates
		### this AS is NOT conserved AS event
                if (not defined $conserved_AS{$Gene_ID}{$A_start}{$A_end}) {
			if (not defined $output{$Chr}{$A_start}{$A_end}) {
                                print "$exon_span_size\t$species_name\tNo\n";
                                $output{$Chr}{$A_start}{$A_end}++;
                        }
		}
		### this AS is YES conserved AS event
                if (defined $conserved_AS{$Gene_ID}{$A_start}{$A_end}) {
	 		if (not defined $output{$Chr}{$A_start}{$A_end}) {
				print "$exon_span_size\t$species_name\tYes\n";
				$output{$Chr}{$A_start}{$A_end}++;
			}
		}
}
close AS;
