#!/usr/bin/perl -w
use strict;
### this script calculate the number of events and genes from the clusters for each type in the second column for "ten species" seperately after remove w22
### ClusterID	Num_Of_Species	B73	Sorghum	Rice	Brachypodium	Millet	Banana	Palm	Arabidopsis	Amborella
my $usage = "$0 <conserved_IntronR_10species_20Jan2016_table.txt> <type>";
die $usage unless $#ARGV == 1;
my ($file, $type) = @ARGV;
my (%B73_genes, %Sorghum_genes, %Rice_genes, %Brachypodium_genes, %Millet_genes, %Banana_genes, %Palm_genes, %Arabidopsis_genes, %Amborella_genes);
my (%B73_events, %Sorghum_events, %Rice_events, %Brachypodium_events, %Millet_events, %Banana_events, %Palm_events, %Arabidopsis_events, %Amborella_events);
my $B73_events_count = 0;
my $Sorghum_events_count = 0;
my $Rice_events_count = 0;
my $Brachypodium_events_count = 0;
my $Millet_events_count = 0;
my $Banana_events_count = 0;
my $Palm_events_count = 0;
my $Arabidopsis_events_count = 0;
my $Amborella_events_count = 0;

my $B73_genes_count = 0;
my $Sorghum_genes_count = 0;
my $Rice_genes_count = 0;
my $Brachypodium_genes_count = 0;
my $Millet_genes_count = 0;
my $Banana_genes_count = 0;
my $Palm_genes_count = 0;
my $Arabidopsis_genes_count = 0;
my $Amborella_genes_count = 0;

my $B73_genelist_output = "B73_genelist_".$type."_".$file;
my $Sorghum_genelist_output = "Sorghum_genelist_".$type."_".$file;
my $Rice_genelist_output = "Rice_genelist_".$type."_".$file;
my $Brachypodium_genelist_output = "Brachypodium_genelist_".$type."_".$file;
my $Millet_genelist_output = "Millet_genelist_".$type."_".$file;
my $Banana_genelist_output = "Banana_genelist_".$type."_".$file;
my $Palm_genelist_output = "Palm_genelist_".$type."_".$file;
my $Arabidopsis_genelist_output = "Arabidopsis_genelist_".$type."_".$file;
my $Amborella_genelist_output = "Amborella_genelist_".$type."_".$file;

my $B73_event_output = "B73_event_".$type."_".$file;
my $Sorghum_event_output = "Sorghum_event_".$type."_".$file;
my $Rice_event_output = "Rice_event_".$type."_".$file;
my $Brachypodium_event_output = "Brachypodium_event_".$type."_".$file;
my $Millet_event_output = "Millet_event_".$type."_".$file;
my $Banana_event_output = "Banana_event_".$type."_".$file;
my $Palm_event_output = "Palm_event_".$type."_".$file;
my $Arabidopsis_event_output = "Arabidopsis_event_".$type."_".$file;
my $Amborella_event_output = "Amborella_event_".$type."_".$file;

open (B73_genelist, ">$B73_genelist_output")  || die "Can't open $B73_genelist_output: $!\n";
open (Sorghum_genelist, ">$Sorghum_genelist_output")  || die "Can't open $Sorghum_genelist_output: $!\n";
open (Rice_genelist, ">$Rice_genelist_output")  || die "Can't open $Rice_genelist_output: $!\n";
open (Brachypodium_genelist, ">$Brachypodium_genelist_output")  || die "Can't open $Brachypodium_genelist_output: $!\n";
open (Millet_genelist, ">$Millet_genelist_output")  || die "Can't open $Millet_genelist_output: $!\n";
open (Banana_genelist, ">$Banana_genelist_output")  || die "Can't open $Banana_genelist_output: $!\n";
open (Palm_genelist, ">$Palm_genelist_output")  || die "Can't open $Palm_genelist_output: $!\n";
open (Arabidopsis_genelist, ">$Arabidopsis_genelist_output")  || die "Can't open $Arabidopsis_genelist_output: $!\n";
open (Amborella_genelist, ">$Amborella_genelist_output")  || die "Can't open $Amborella_genelist_output: $!\n";

open (B73_event, ">$B73_event_output")  || die "Can't open $B73_event_output: $!\n";
open (Sorghum_event, ">$Sorghum_event_output")  || die "Can't open $Sorghum_event_output: $!\n";
open (Rice_event, ">$Rice_event_output")  || die "Can't open $Rice_event_output: $!\n";
open (Brachypodium_event, ">$Brachypodium_event_output")  || die "Can't open $Brachypodium_event_output: $!\n";
open (Millet_event, ">$Millet_event_output")  || die "Can't open $Millet_event_output: $!\n";
open (Banana_event, ">$Banana_event_output")  || die "Can't open $Banana_event_output: $!\n";
open (Palm_event, ">$Palm_event_output")  || die "Can't open $Palm_event_output: $!\n";
open (Arabidopsis_event, ">$Arabidopsis_event_output")  || die "Can't open $Arabidopsis_event_output: $!\n";
open (Amborella_event, ">$Amborella_event_output")  || die "Can't open $Amborella_event_output: $!\n";

open (FILE, $file) or die "Can not open input file: $file\n";
while (<FILE>) {
			chomp;
			next if /^ASeventID/; 		
			my ($ClusterID, $Num_Of_Species, $B73, $Sorghum, $Rice, $Brachypodium, $Millet, $Banana, $Palm, $Arabidopsis, $Amborella) = split (/\s+/, $_);
			if ($_ =~ /$type/) {
				my @array_B73 = split (/\;/, $B73);
				foreach my $i (@array_B73) {
					next if $B73 =~ /^NA/;
					my ($B73_gene, $location) = split (/\:/, $i);
					if (not defined $B73_events{$i} && $i !~ /^NA/){
							$B73_events_count++;
						    	$B73_events{$i}++;
							print B73_event "$i\n";
					}
					if (not defined $B73_genes{$B73_gene} && $i !~ /^NA/){
							$B73_genes_count++;
						    	$B73_genes{$B73_gene}++;
						    	print B73_genelist "$B73_gene\n";
					}
				}
				my @array_Sorghum = split (/\;/, $Sorghum);
                                foreach my $i (@array_Sorghum) {
					next if $Sorghum =~ /^NA/;
                                        my ($Sorghum_gene, $location) = split (/\:/, $i);
                                        if (not defined $Sorghum_events{$i} && $i !~ /^NA/){
                                                        $Sorghum_events_count++;
                                                        $Sorghum_events{$i}++;
							print Sorghum_event "$i\n";
                                        }
                                        if (not defined $Sorghum_genes{$Sorghum_gene} && $i !~ /^NA/){
                                                        $Sorghum_genes_count++;
                                                        $Sorghum_genes{$Sorghum_gene}++;
                                                        print Sorghum_genelist "$Sorghum_gene\n";
                                        }
                                }
				my @array_Rice = split (/\;/, $Rice);
                                foreach my $i (@array_Rice) {
					next if $Rice =~ /^NA/;
                                        my ($Rice_gene, $location) = split (/\:/, $i);
                                        if (not defined $Rice_events{$i} && $i !~ /^NA/){
                                                        $Rice_events_count++;
                                                        $Rice_events{$i}++;
							print Rice_event "$i\n";
                                        }
                                        if (not defined $Rice_genes{$Rice_gene} && $i !~ /^NA/){
                                                        $Rice_genes_count++;
                                                        $Rice_genes{$Rice_gene}++;
                                                        print Rice_genelist "$Rice_gene\n";
                                        }
                                }
				my @array_Brachypodium = split (/\;/, $Brachypodium);
                                foreach my $i (@array_Brachypodium) {
					next if $Brachypodium =~ /^NA/;
                                        my ($Brachypodium_gene, $location) = split (/\:/, $i);
                                        if (not defined $Brachypodium_events{$i} && $i !~ /^NA/){
                                                        $Brachypodium_events_count++;
                                                        $Brachypodium_events{$i}++;
							print Brachypodium_event "$i\n";
                                        }
                                        if (not defined $Brachypodium_genes{$Brachypodium_gene} && $i !~ /^NA/){
                                                        $Brachypodium_genes_count++;
                                                        $Brachypodium_genes{$Brachypodium_gene}++;
                                                        print Brachypodium_genelist "$Brachypodium_gene\n";
                                        }
                                }
				my @array_Millet = split (/\;/, $Millet);
                                foreach my $i (@array_Millet) {
					next if $Millet =~ /^NA/;
                                        my ($Millet_gene, $location) = split (/\:/, $i);
                                        if (not defined $Millet_events{$i} && $i !~ /^NA/){
                                                        $Millet_events_count++;
                                                        $Millet_events{$i}++;
							print Millet_event "$i\n";
                                        }
                                        if (not defined $Millet_genes{$Millet_gene} && $i !~ /^NA/){
                                                        $Millet_genes_count++;
                                                        $Millet_genes{$Millet_gene}++;
                                                        print Millet_genelist "$Millet_gene\n";
                                        }
                                }
				my @array_Banana = split (/\;/, $Banana);
                                foreach my $i (@array_Banana) {
					next if $Banana =~ /^NA/;
                                        my ($Banana_gene, $location) = split (/\:/, $i);
                                        if (not defined $Banana_events{$i} && $i !~ /^NA/){
                                                        $Banana_events_count++;
                                                        $Banana_events{$i}++;
							print Banana_event "$i\n";
                                        }
                                        if (not defined $Banana_genes{$Banana_gene} && $i !~ /^NA/){
                                                        $Banana_genes_count++;
                                                        $Banana_genes{$Banana_gene}++;
                                                        print Banana_genelist "$Banana_gene\n";
                                        }
                                }
				my @array_Palm = split (/\;/, $Palm);
                                foreach my $i (@array_Palm) {
					next if $Palm =~ /^NA/;
                                        my ($Palm_gene, $location) = split (/\:/, $i);
                                        if (not defined $Palm_events{$i} && $i !~ /^NA/){
                                                        $Palm_events_count++;
                                                        $Palm_events{$i}++;
							print Palm_event "$i\n";
                                        }
                                        if (not defined $Palm_genes{$Palm_gene} && $i !~ /^NA/){
                                                        $Palm_genes_count++;
                                                        $Palm_genes{$Palm_gene}++;
                                                        print Palm_genelist "$Palm_gene\n";
                                        }
                                }
				my @array_Arabidopsis = split (/\;/, $Arabidopsis);
                                foreach my $i (@array_Arabidopsis) {
					next if $Arabidopsis =~ /^NA/;
                                        my ($Arabidopsis_gene, $location) = split (/\:/, $i);
                                        if (not defined $Arabidopsis_events{$i} && $i !~ /^NA/){
                                                        $Arabidopsis_events_count++;
                                                        $Arabidopsis_events{$i}++;
							print Arabidopsis_event "$i\n";
                                        }
                                        if (not defined $Arabidopsis_genes{$Arabidopsis_gene} && $i !~ /^NA/){
                                                        $Arabidopsis_genes_count++;
                                                        $Arabidopsis_genes{$Arabidopsis_gene}++;
                                                        print Arabidopsis_genelist "$Arabidopsis_gene\n";
                                        }
                                }
				my @array_Amborella = split (/\;/, $Amborella);
                                foreach my $i (@array_Amborella) {
					next if $Amborella =~ /^NA/;
                                        my ($Amborella_gene, $location) = split (/\:/, $i);
                                        if (not defined $Amborella_events{$i}){
                                                        $Amborella_events_count++;
                                                        $Amborella_events{$i}++;
							print Amborella_event "$i\n";
                                        }
                                        if (not defined $Amborella_genes{$Amborella_gene}){
                                                        $Amborella_genes_count++;
                                                        $Amborella_genes{$Amborella_gene}++;
                                                        print Amborella_genelist "$Amborella_gene\n";
                                        }
                                }	
			}
}
close FILE;

print "B73	events_count:$B73_events_count	genes_count:$B73_genes_count\n";
print "Sorghum	events_count:$Sorghum_events_count	genes_count:$Sorghum_genes_count\n";
print "Rice	events_count:$Rice_events_count	genes_count:$Rice_genes_count\n";
print "Brachypodium	events_count:$Brachypodium_events_count	genes_count:$Brachypodium_genes_count\n";
print "Millet	events_count:$Millet_events_count	genes_count:$Millet_genes_count\n";
print "Banana	events_count:$Banana_events_count	genes_count:$Banana_genes_count\n";
print "Palm	events_count:$Palm_events_count	genes_count:$Palm_genes_count\n";
print "Arabidopsis	events_count:$Arabidopsis_events_count	genes_count:$Arabidopsis_genes_count\n";
print "Amborella	events_count:$Amborella_events_count	genes_count:$Amborella_genes_count\n";






            
