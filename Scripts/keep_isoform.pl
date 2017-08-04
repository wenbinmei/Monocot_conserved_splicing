#!/usr/bin/perl -w
use strict;

### keep only the isoform ID after filtering FPKM and isoform ratio

my $usage = "$0 <Feb6_Mo17_CombinedData_coverage2_Frac0.05_ASloci_overlapMaizeGene.gtf> <B73_isoID_afterFilter_FPKMIsoRatio>";

die $usage unless $#ARGV == 1;

#chr1    assembler       transcript      4686    9651    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";
#chr1    assembler       exon    4686    5188    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";
#chr1    assembler       exon    5342    5407    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";
#chr1    assembler       exon    5857    5975    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";
#chr1    assembler       exon    6108    6265    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";
#chr1    assembler       exon    6362    6517    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";
#chr1    assembler       exon    6639    6797    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";
#chr1    assembler       exon    6918    7120    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";
#chr1    assembler       exon    7594    7903    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";
#chr1    assembler       exon    9193    9651    .       -       .       gene_id "PASA_cluster_1"; transcript_id "align_id:1142719|asmbl_1";

my ($GTF, $FILE) = @ARGV;
my %keep_IR_iso;

open (File, $FILE) or die "Can not open input file: $FILE\n";
while (<File>) {
          chomp;
          if (not defined $keep_IR_iso{$_}) {
                    $keep_IR_iso{$_}++;
          }
}
close File;


open (Transcript, $GTF) or die "Can not open input file: $GTF\n";
while (<Transcript>) {
         chomp;
         next if /^\s*$/;
         my $info = $_;
         my @line = split (/\s+/, $_);
         my $transID = $line[11];
         $transID =~ s/[";]//g;
         #print "$transID\n"; 
         if (defined $keep_IR_iso{$transID}) {
                 print "$_\n";
         }
}
close Transcript;
