#!/usr/bin/perl -w
use strict;

###################################################################
# Author: Wenbin Mei
# 
# This scripts take a filtered GFF file and identify the overlap gene name from annotation, and only keep the isoforms with five different alternative splicing events;
# filter the PASA_OUTPUT_pasa_wmei_B73_March18th_2015.indiv_splice_labels_and_coords.dat file and only keep the events both A_isoform and B_isoform exist in the GTF file;
#
# step1: identify the gene in the cuffcompare file and store the isoform, gene information and make sure each gene have at least two isoforms 
# 
# step2: process the PASA output file and store the AS isoform list
#
# step3: store the final sets of isoform ID only in those genes with at least two isoforms
#        here, the difference between genes sets also with two isoforms in step1 is, in step1 you do not know whether these two isoforms are involves in five types of AS or not, here we do !!!
#
# step4: filter the PASA output file and only keep those isoforms ID in those genes with at least two isoforms from step3
#
# step5: filter the final GTF file, you can see the num of isoforms in final PASA OUTPUT is the same as the final GTF file
#
#
##################################################################



# File1 (Filtered GTF file)
# chr1    assembler       transcript      62297   63517   .       -       .       gene_id "PASA_cluster_3"; transcript_id "align_id:1142739|asmbl_21";
# chr1    assembler       exon    62297   62371   .       -       .       gene_id "PASA_cluster_3"; transcript_id "align_id:1142739|asmbl_21";
# chr1    assembler       exon    62542   63517   .       -       .       gene_id "PASA_cluster_3"; transcript_id "align_id:1142739|asmbl_21";
# chr1    assembler       transcript      71772   72061   .       -       .       gene_id "PASA_cluster_3"; transcript_id "align_id:1142741|asmbl_23";
# chr1    assembler       exon    71772   71951   .       -       .       gene_id "PASA_cluster_3"; transcript_id "align_id:1142741|asmbl_23";
# chr1    assembler       exon    72033   72061   .       -       .       gene_id "PASA_cluster_3"; transcript_id "align_id:1142741|asmbl_23";
# chr1    assembler       transcript      84675   93563   .       -       .       gene_id "PASA_cluster_4"; transcript_id "align_id:1142745|asmbl_27";
# chr1    assembler       exon    84675   85000   .       -       .       gene_id "PASA_cluster_4"; transcript_id "align_id:1142745|asmbl_27";
# chr1    assembler       exon    86014   86169   .       -       .       gene_id "PASA_cluster_4"; transcript_id "align_id:1142745|asmbl_27";
# chr1    assembler       exon    90685   90798   .       -       .       gene_id "PASA_cluster_4"; transcript_id "align_id:1142745|asmbl_27";
# chr1    assembler       exon    92521   92622   .       -       .       gene_id "PASA_cluster_4"; transcript_id "align_id:1142745|asmbl_27";
# chr1    assembler       exon    93386   93563   .       -       .       gene_id "PASA_cluster_4"; transcript_id "align_id:1142745|asmbl_27";


# File2 (cuffcompare file with maize geneID)
# ref_gene_id     ref_id  class_code      cuff_gene_id    cuff_id FMI     FPKM    FPKM_conf_lo    FPKM_conf_hi    cov     len     major_iso_id    ref_match_len
# -       -       u       S1      asmbl_1 0       0.000000        0.000000        0.000000        0.000000        3165    asmbl_3 -
# -       -       u       S1      asmbl_4 0       0.000000        0.000000        0.000000        0.000000        2754    asmbl_3 -
# -       -       u       S1      asmbl_3 0       0.000000        0.000000        0.000000        0.000000        3552    asmbl_3 -
# -       -       u       S1      asmbl_2 0       0.000000        0.000000        0.000000        0.000000        2646    asmbl_3 -
# -       -       u       S1      asmbl_6 0       0.000000        0.000000        0.000000        0.000000        2728    asmbl_3 -
# -       -       u       S1      asmbl_5 0       0.000000        0.000000        0.000000        0.000000        2858    asmbl_3 -
# GRMZM2G059865   GRMZM2G059865_T02       c       S2      asmbl_8 0       0.000000        0.000000        0.000000        0.000000        1968    asmbl_8 2412
# GRMZM2G059865   GRMZM2G059865_T01       =       S2      asmbl_7 0       0.000000        0.000000        0.000000        0.000000        1938    asmbl_8 1966


# File3 PASA_OUTPUT_pasa_wmei_B73_March18th_2015.indiv_splice_labels_and_coords.dat
# Chr     AS_type strand  Gene_ID Isoform_A       A_start A_end   Isoform_B       B_start B_end
# chr1    skipped_exon    -       1       asmbl_1 5857    5407    asmbl_3 5694    5515
# chr1    alt_donor       -       1       asmbl_1 5857    5407    asmbl_4 5515    5407
# chr1    retained_intron -       1       asmbl_2 NA      NA      asmbl_3 5857    5694
# chr1    retained_intron -       1       asmbl_4 NA      NA      asmbl_3 5857    5694
# chr1    retained_intron -       2       asmbl_10        NA      NA      asmbl_11        48952   48564
# chr1    retained_intron -       2       asmbl_11        NA      NA      asmbl_10        47636   46679


my $usage = "$0 <Filtered GTF File> <cuffcompare file> <PASA OUTPUT file>";
die $usage unless $#ARGV == 2;
my (%gene_list, %GN_two_transcript, %AS_isoform);
my $GN = " ";
#my ($count, $j);
my ($Filtered_GTF, $asmbl_annotation, $PASA_OUTPUT) = @ARGV;

### step1
### open the file link the PASA isoform with maize annotation, make sure they still have more than two isoform per gene
my %count_transcript;
open (CUFFCOMPARE, $asmbl_annotation) or die "Can not open input file: $asmbl_annotation\n";
while (<CUFFCOMPARE>) {
            chomp;
             next if /ref_gene_id/;   ### skip the header line
             my @line = split (/\t/, $_);
             $GN = $line[0];          ### gene name
             next if ($GN =~ /-/);    ### no annotation association gene name in the cuffcompare
             my $Iso_id = $line[4];
             my @list = split (/\|/, $Iso_id);   ### parse align_id:1142745|asmbl_27
                $Iso_id = $list[1];
             $count_transcript{$GN}++;
             $gene_list{$Iso_id} = $GN;
             if ( (not defined $GN_two_transcript{$GN}) && ($count_transcript{$GN} > 1) ) { ### make sure each gene still have at least two isoforms
                                $GN_two_transcript{$GN} = 0;
             }
}
close CUFFCOMPARE;

### step2
open (PASA, $PASA_OUTPUT) or die "Can not open input file: $PASA_OUTPUT\n";
while (<PASA>) {
             chomp;
             if ($_ =~ /AS_type/) {
                           print "$_\n";
                           next;
             }
             my @line = split (/\s+/, $_);
             my $A_isoform = $line[4];
             my $B_isoform = $line[7];
             if (not defined $AS_isoform{$A_isoform}) {
                           $AS_isoform{$A_isoform}++;
             }
             if (not defined $AS_isoform{$B_isoform}) {
                           $AS_isoform{$B_isoform}++;
             }
}
close PASA;       

my %final_iso;

### step3
open (FILE, $Filtered_GTF) or die "Can not open input file: $Filtered_GTF\n";
while (<FILE>) {
            chomp;
            if (!/^\s*$/) {   ### remove the blank line
               my @line = split (/\s+/, $_);
               my $type = $line[2];
               if ( $_ =~ /gene_id\s+"(\S+)";\s+transcript_id\s+"(\S+)";/) {
                    my $gene = $1;
                    my $transID = $2;
                    my @list = split (/\|/, $transID);   ### parse align_id:1142745|asmbl_27
                    my $asmbl_ID = $list[1];
                    if ( (defined $gene_list{$asmbl_ID}) && (defined $AS_isoform{$asmbl_ID}) ) {
                              my $MaizeGene = $gene_list{$asmbl_ID};
                              if (defined $GN_two_transcript{$MaizeGene}) {
                                   my $line = $_;
                                   $final_iso{$asmbl_ID} = $MaizeGene;
                                   $asmbl_ID = $MaizeGene."_".$asmbl_ID;
                                   $line =~ s/$gene/$MaizeGene/;
                                   $line =~ s/\|/\_/;               ### original $transID has the "|" which means "or", so I first replace "|" with "_";
                                   $transID =~ s/\|/\_/;
                                   $line =~ s/$transID/$asmbl_ID/;
                                   #print "$transID\n";
                                   #print "$asmbl_ID\n";
                                   print "$line\n";
                              }
                    }
               }
            }
}
close FILE;

#my $cmd  = 'perl /home/wmei/lfs/PROJECT/FROM_SCRATCH/src/Alternative_Splicing_Pipeline_Scripts/filter_single_transcript_loci_PASA.pl';
#   $cmd .= ' temp.gtf';
#   $cmd .= ' > test.gtf';
#system($cmd);

my %gtf_keep_iso;

open(OUTPUT, ">Final_$PASA_OUTPUT");

### step4
open (PASA_OUTPUT, $PASA_OUTPUT) or die "Can not open input file: $PASA_OUTPUT\n";
while (<PASA_OUTPUT>) {
             chomp;
             print "$_\n";
             if ($_ =~ /AS_type/) {
                    print OUTPUT "$_\n";
                    next;
             }
             my @line = split (/\s+/, $_);
             #print OUTPUT "$line[4]\n";
             #print OUTPUT "$line[7]\n";
             if ( (defined $final_iso{$line[4]}) && (defined $final_iso{$line[7]}) && ( $final_iso{$line[4]} eq $final_iso{$line[7]}) ) {
                    #Chr	AS_type	strand	Gene_ID	Isoform_A	A_start	A_end	Isoform_B	B_start	B_end
                    my $output_gene = $line[3];
                    my $A_isoform = $line[4];
                    my $B_isoform = $line[7];
                    #print "$A_isoform\n";
                    #print "$gene_list{$A_isoform} /n";
                    $output_gene =~ s/$line[3]/$gene_list{$A_isoform}/;
                    print OUTPUT "$line[0]\t$line[1]\t$line[2]\t$output_gene\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\n";
                    $gtf_keep_iso{$A_isoform}++;
                    $gtf_keep_iso{$B_isoform}++;
             }
             elsif (not defined $final_iso{$line[4]}) {
                    next;
             }
             elsif (not defined $final_iso{$line[7]}) {
                    next;
             }
}
close PASA_OUTPUT;

open(OUTPUT_GTF, ">Final_$Filtered_GTF");

### step5
open (FILE_OUTPUT, $Filtered_GTF) or die "Can not open input file: $Filtered_GTF\n";
while (<FILE_OUTPUT>) {
            chomp;
            if (!/^\s*$/) {
            	my @line = split (/\s+/, $_);
            	my $type = $line[2];
            	if ( $_ =~ /gene_id\s+"(\S+)";\s+transcript_id\s+"(\S+)";/) {
                	my $gene = $1;
                    	my $transID = $2;
                    	my @list = split (/\|/, $transID);
                    	my $asmbl_ID = $list[1];
                        my $MaizeGene = $gene_list{$asmbl_ID};
                        if (defined $gtf_keep_iso{$asmbl_ID}) {
                            my $line = $_;
                            $asmbl_ID = $MaizeGene."_".$asmbl_ID;
                            $line =~ s/$gene/$MaizeGene/;    
                            $line =~ s/\|/\_/;               ### original $transID has the "|" which means "or", so I first replace "|" with "_";
                            $transID =~ s/\|/\_/;
                            $line =~ s/$transID/$asmbl_ID/;
                            print OUTPUT_GTF "$line\n";
                        }

            	}
            }
}
close FILE_OUTPUT;

#system('rm temp.gtf');
#system('rm test.gtf');




































