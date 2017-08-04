#!/usr/bin/perl -w
use strict;

### this script take a list of genes and the AS output file generate previous and identify the subsets the AS events match the gene list

my $usage = "$0 <AS event> <gene list>";
die $usage unless $#ARGV == 1;
my ($AS_event, $gene_list) = @ARGV;
my %keep;

open (GENE, $gene_list) or die "Can not open input file: $gene_list\n";
while (<GENE>) {
          chomp;
          my @line = split(/\s+/, $_);
          $keep{$line[0]} = $line[0];
}
close GENE;

open (AS, $AS_event) or die "Can not open input file: $AS_event\n";
while (<AS>) {
          chomp;
          my ($chr, $AS_type, $strand, $Gene_ID, $Isoform_A, $A_start, $A_end, $Isoform_B, $B_start, $B_end) = split(/\s+/, $_);
          if (defined $keep{$Gene_ID}) {
		print "$_\n";
	  }
}
close AS;

