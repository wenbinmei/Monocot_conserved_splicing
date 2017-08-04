#!/usr/bin/perl -w
use strict;

### filter out the PASA GTF file already filtered Isoform Fraction 0.05 and 2 reads support junction early, Now filtering additional IR isoform with low coverage

my $usage = "$0 <Feb6_Mo17_CombinedData_coverage2_Frac0.05_ASloci_overlapMaizeGene.gtf> <low coverage IR isoform>";

die $usage unless $#ARGV == 1;

#chr1	assembler	transcript	109505	111136	.	-	.	gene_id "GRMZM2G093344"; transcript_id "align_id:706660|asmbl_78";
#chr1	assembler	transcript	108940	110529	.	-	.	gene_id "GRMZM2G093344"; transcript_id "align_id:706661|asmbl_79";
#chr1	assembler	transcript	143716	145635	.	-	.	gene_id "GRMZM5G833153"; transcript_id "align_id:706662|asmbl_80";
#chr1	assembler	transcript	143716	145609	.	-	.	gene_id "GRMZM5G833153"; transcript_id "align_id:706664|asmbl_82";
#chr1	assembler	transcript	515984	519683	.	-	.	gene_id "GRMZM2G330436"; transcript_id "align_id:706679|asmbl_97";
#chr1	assembler	transcript	515997	522398	.	-	.	gene_id "GRMZM2G330436"; transcript_id "align_id:706680|asmbl_98";
#chr1	assembler	transcript	515997	522398	.	-	.	gene_id "GRMZM2G330436"; transcript_id "align_id:706681|asmbl_99";
#chr1	assembler	transcript	644477	658856	.	-	.	gene_id "GRMZM2G032104"; transcript_id "align_id:706685|asmbl_103";
#chr1	assembler	transcript	644477	650249	.	-	.	gene_id "GRMZM2G032104"; transcript_id "align_id:706689|asmbl_107";
#chr1	assembler	transcript	644477	658828	.	-	.	gene_id "GRMZM2G032104"; transcript_id "align_id:706701|asmbl_119";

my ($GTF, $FILE) = @ARGV;
my %filter_IR_iso;

open (File, $FILE) or die "Can not open input file: $FILE\n";
while (<File>) {
          chomp;
          if (not defined $filter_IR_iso{$_}) {
                    $filter_IR_iso{$_}++;    ### a list of low coverage IR isoform
          }
}
close File;

#### The list of isoforms in the GTF file has been filtered out using FPKM >= 1 and 2 reads support the junction, just need to merge the IR isoform after filtering

open (Transcript, $GTF) or die "Can not open input file: $GTF\n";
while (<Transcript>) {
         chomp;
         next if /^\s*$/;
         my $info = $_;
         my @line = split (/\s+/, $_);
         my $transID = $line[11];  
         if ($transID =~ /(asmbl_\d+)/) {
              $transID = $1;
         }
              
         if (not defined $filter_IR_iso{$transID}) {
                 print "$_\n";
         }
}
close Transcript;
