#!/usr/bin/perl -w
use strict;
my $usage = "$0 <conserved_splicing_file>";
die $usage unless $#ARGV == 0;
my ($file) = @ARGV;

### this script identified the conserved clustered:
# 1) at least two of the three species in PACMAD conserved (Maize, Sorghum, Millet)
# 2) at least one species in PACMAD and one species in BEP (Rice, Brachypodium) conserved
### ClusterID   Num_Of_Species  B73	Sorghum Rice    Brachypodium    Millet  Banana  Palm    Arabidopsis     Amborella

my (%B73_genes, %Sorghum_genes, %Rice_genes, %Brachypodium_genes, %Millet_genes, %Banana_genes, %Palm_genes, %Arabidopsis_genes, %Amborella_genes);
my $conserved_PACMAD = "conserved_ASevents_PACMAD.csv";
my $conserved_BEP = "conserved_ASevents_BEP.csv";
my $conserved_BEP_PACMAD = "conserved_ASevents_conserved_BEP_PACMAD.csv";
my $all_conserved_grass = "all_conserved_grass.csv";

open (PACMAD, ">$conserved_PACMAD")  || die "Can't open $conserved_PACMAD: $!\n";
open (BEP, ">$conserved_BEP")  || die "Can't open $conserved_BEP: $!\n";
open (BEP_PACMAD, ">$conserved_BEP_PACMAD")  || die "Can't open $conserved_BEP_PACMAD: $!\n";
open (ALL_CONSERVED, ">$all_conserved_grass")  || die "Can't open $all_conserved_grass: $!\n";

open (FILE, $file) or die "Can not open input file: $file\n";
while (<FILE>) {
                        chomp;
			my $B73_count = 0;
			my $Maize_count = 0;
			my $Sorghum_count = 0;
			my $Rice_count	= 0;
			my $Brachypodium_count = 0;
			my $Millet_count = 0;
                        next if /^ASeventID/;           
                        my ($ClusterID, $Num_Of_Species, $B73, $Sorghum, $Rice, $Brachypodium, $Millet, $Banana, $Palm, $Arabidopsis, $Amborella) = split (/\s+/, $_);
			if ($B73 !~ /^NA/) {
				$Maize_count = 1;
			}
			if ($Sorghum !~ /^NA/) {
                                $Sorghum_count = 1;
                        }
			if ($Rice !~ /^NA/) {
                                $Rice_count = 1;
                        }
			if ($Brachypodium !~ /^NA/) {
                                $Brachypodium_count = 1;
                        }
			if ($Millet !~ /^NA/) {
                                $Millet_count = 1;
                        }
			my $PACMAD_count = 0;
			my $BEP_count	= 0;
			my $All_count = 0;
			$PACMAD_count = $Maize_count + $Sorghum_count + $Millet_count;
			$BEP_count = $Rice_count + $Brachypodium_count;
			$All_count = $Maize_count + $Sorghum_count + $Millet_count + $Rice_count + $Brachypodium_count;
			### identified conserved in PACMAD (maize, sorghum, millet)
			if ($PACMAD_count >= 2) {
				print PACMAD "$_\n";
			}
			### identified conserved in BEP (Rice, brachypodium)
			if ($BEP_count >= 2) {
                                print BEP "$_\n";
                        }
			### conserved at least one species in PACMAD and ones species in BEP, suggest pan-grass conserved
			if ($PACMAD_count >= 1 && $BEP_count >= 1) {
				print BEP_PACMAD "$_\n";
			}
			### all conserved in grass (Maize, sorghum, rice, foxtail, brachypodium)
			if ($All_count == 5) {
                                print ALL_CONSERVED "$_\n";
                        }
}
close FILE;
