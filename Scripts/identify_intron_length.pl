### this script take two files PASA_AS_File and Conserved_Event_File plus the species name
### and identify the intron length for Intron Retained events which exist in the conserved AS or NOT

#Final_PASA_OUTPUT_pasa_wmei_amborella_Aug2_2015.indiv_splice_labels_and_coords_rmEVM.dat
#AmTr_v1.0_scaffold00001 retained_intron -       AmTr_v1.0_scaffold00001.5       asmbl_13        NA      NA      asmbl_12        82809   81576
#AmTr_v1.0_scaffold00001 skipped_exon    -       AmTr_v1.0_scaffold00001.5       asmbl_12        69194   65019   asmbl_13        65377   65269

# Amborella_event_IntronR_conserved_fiveASevents_10species_20Jan2016_table_RemoveUniqueConservedB73W22.txt
# AmTr_v1.0_scaffold00096.43:IntronR_AmTr_v1.0_scaffold00096_+_NA_NA_1064072_1064258
# AmTr_v1.0_scaffold00024.131:IntronR_AmTr_v1.0_scaffold00024_-_NA_NA_2821012_2819791

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
		my $start = $line[-2];
		my $end = $line[-1];
		#print "$chr\t$strand\t$start\t$end\n";
		if (not defined $conserved_AS{$gene}{$start}{$end}) {
			$conserved_AS{$gene}{$start}{$end} = $gene."_".$start."_".$end;
		}
}
close CONSERVED;

print "Intron_length\tSpecies\tConserved\n";

open (AS, $AS_file) or die "Can not open input file: $AS_file\n";
while (<AS>) {
		next if ($_ !~ /retained_intron/);
		my ($Chr, $AS_type, $strand, $Gene_ID, $Isoform_A, $A_start, $A_end, $Isoform_B, $B_start, $B_end) = split (/\s+/, $_);
                my $intron_size = abs($B_start - $B_end) - 1; ### because the coordinate is based on exon coordinates
		### this AS is NOT conserved AS event
                if (not defined $conserved_AS{$Gene_ID}{$B_start}{$B_end}) {
			if (not defined $output{$Chr}{$B_start}{$B_end}) {
                                print "$intron_size\t$species_name\tNo\n";
                                $output{$Chr}{$B_start}{$B_end}++;
                        }
		}
		### this AS is YES conserved AS event
                if (defined $conserved_AS{$Gene_ID}{$B_start}{$B_end}) {
	 		if (not defined $output{$Chr}{$B_start}{$B_end}) {
				print "$intron_size\t$species_name\tYes\n";
				$output{$Chr}{$B_start}{$B_end}++;
			}
		}
}
close AS;
