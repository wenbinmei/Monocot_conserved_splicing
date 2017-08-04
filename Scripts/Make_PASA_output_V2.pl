#!/usr/bin/perl -w
use strict;
################################
#  Feb10/2014 author: Wenbin Mei  I change the scripts that the output is skipped exon not retained exon.
#  Aug4/2015  author: Wenbin Mei  additonal change based on V1, 1) all the coordinate is based on the exon not intron, so for example the coordinate of intron retention is the first bp of two adjacent exon;
#                                 2) have four different coordinate for a pair of alternative splicing events, details are illustrated below
#                                 3) always the longer distance is isoformA relative to isoformB, it would be easier for downstream to pick up the fasta junction seq                          
#
#   =: exon; -: intron; #: retained part in the isoforms
#              3       4                                    3      4       
#    IR: =======-------==========        ES: =======--------========-------========
#        =======#######==========            =======-----------------------========                                                                                                          
#              1       2                           1                       2
#
#
#                   3        4                              3           4
#    AltD:  =========--------==========      AltA: ==========-----------==============
#           ======-----------==========            ==========---------------==========
#                1           2                              1               2
#
#                                              3        4
#    AtTE:  =====---------=======-------========--------======
#           =====---------=======-------========--------------------------=========
#                                              1                          2 
#                                              
#################################

### This scripts takes the PASA output file xxxxxx.indiv_splice_labels_and_coords.dat and rewrite in a format which is better to identify the information for the downstream analysis

my $usage = "$0 <pasa.indiv_splice_labels_and_coords.dat> <pasa_assemblies.gtf>";

#AmTr_v1.0_scaffold00001 asmbl_183671    5       alt_acceptor    61845   61846   -       asmbl_31
#AmTr_v1.0_scaffold00001 asmbl_183672    5       alt_acceptor    61915   61916   -       asmbl_36

#AmTr_v1.0_scaffold00001 asmbl_183736    46      alt_donor       334793  334794  +       asmbl_146
#AmTr_v1.0_scaffold00001 asmbl_183740    46      alt_donor       335324  335325  +       asmbl_168

#AmTr_v1.0_scaffold00001 asmbl_183673    5       retained_intron 81577   82636   -       asmbl_10
#AmTr_v1.0_scaffold00001 asmbl_183674    5       spliced_intron  81577   82636   -       asmbl_11

#AmTr_v1.0_scaffold00001 asmbl_183681    5       retained_exon   65269   65377   -       asmbl_9
#AmTr_v1.0_scaffold00001 asmbl_183673    5       skipped_exon    65020   69193   -       asmbl_10

#Chr1    asmbl_15003     14997   alternate_exon  62159   62216   +       TCONS_00000015
#Chr1    asmbl_15004     14997   alternate_exon  62050   62054   +       TCONS_00000018

### alt_acceptor, alt_donor and alternate_exon two isoforms are next to each other
### retained_intron and spliced_intron are a pair
### retained_exon and skipped_exon are a pair

my ($File, $PASA_GTF) = @ARGV;
my $i = 0;
my (%link, %list_transID, %exon_start, %exon_end);

### only five different AS types, and filter the rest of other events and all these events are in pair and two pairs are next to each other;
open (FILE, $File) or die "Can not open input file: $File\n";
while (<FILE>) {
        chomp;
        my ($chr, $isoID, $geneID, $type, $start, $end, $strand, $cDNA_accs) = split (/\s+/, $_);
        if ($type =~ /(alt_acceptor|alt_donor|retained_exon|retained_intron|skipped_exon|spliced_intron|alternate_exon)/) {
               $i++;
               $link{$i} = $_;
        }
}
close FILE;

### open PASA gtf file and store the exon start and end location in a hash

#Chr1    assembler       transcript      3630    6263    .       +       .       gene_id "PASA_cluster_14994"; transcript_id "align_id:873078|asmbl_14994";
#Chr1    assembler       exon    3630    3913    .       +       .       gene_id "PASA_cluster_14994"; transcript_id "align_id:873078|asmbl_14994";
#Chr1    assembler       exon    3996    4276    .       +       .       gene_id "PASA_cluster_14994"; transcript_id "align_id:873078|asmbl_14994";
#Chr1    assembler       exon    4486    4605    .       +       .       gene_id "PASA_cluster_14994"; transcript_id "align_id:873078|asmbl_14994";
#Chr1    assembler       exon    4706    5095    .       +       .       gene_id "PASA_cluster_14994"; transcript_id "align_id:873078|asmbl_14994";
#Chr1    assembler       exon    5174    5326    .       +       .       gene_id "PASA_cluster_14994"; transcript_id "align_id:873078|asmbl_14994";
#Chr1    assembler       exon    5439    6263    .       +       .       gene_id "PASA_cluster_14994"; transcript_id "align_id:873078|asmbl_14994";

#Chr1    assembler       transcript      50012   51118   .       -       .       gene_id "PASA_cluster_17966"; transcript_id "align_id:881583|asmbl_23499";
#Chr1    assembler       exon    50012   50337   .       -       .       gene_id "PASA_cluster_17966"; transcript_id "align_id:881583|asmbl_23499";
#Chr1    assembler       exon    50419   50447   .       -       .       gene_id "PASA_cluster_17966"; transcript_id "align_id:881583|asmbl_23499";
#Chr1    assembler       exon    50496   50631   .       -       .       gene_id "PASA_cluster_17966"; transcript_id "align_id:881583|asmbl_23499";
#Chr1    assembler       exon    50883   50963   .       -       .       gene_id "PASA_cluster_17966"; transcript_id "align_id:881583|asmbl_23499";
#Chr1    assembler       exon    51064   51118   .       -       .       gene_id "PASA_cluster_17966"; transcript_id "align_id:881583|asmbl_23499";

my ($transID, $count);

### store exon location
open (GTF, $PASA_GTF) or die "Can not open input file: $PASA_GTF\n";
while (<GTF>) {
          chomp;
          next if /^\s*$/;
          my @line = split (/\s+/, $_);
          my $type = $line[2];
          my $start = $line[3];
          my $end = $line[4];
          if ( ($type =~ /transcript/) && ($line[11] =~ /(asmbl_\d+)/) ) {
               $transID = $1;
               if (not defined $list_transID{$transID}) {
                          $list_transID{$transID}++;
                          $count = 1;
               }
          }

          if ($type =~ /exon/) {
               $exon_start{$transID}{$count} = $start;
               $exon_end{$transID}{$count}   = $end;
               #print "$exon_start{$transID}{$count}\t$exon_end{$transID}{$count}\t$count\n";
               $count++;
          }
}     
close GTF;                  
                    

open (OUTPUT, '>temp'); 
print OUTPUT "Chr"."\t"."AS_type"."\t"."strand"."\t"."Gene_ID"."\t"."Isoform_A"."\t"."A_start"."\t"."A_end"."\t"."Isoform_B"."\t"."B_start"."\t"."B_end"."\n";

for (my $j = 1; $j <= $i; $j++) {
    if (1 == $j % 2) {   ### this is an odd number
          my ($chr, $isoID, $geneID, $type, $start, $end, $strand, $cDNA_accs) = split (/\s+/, $link{$j});
          print OUTPUT "$chr\t$type\t$strand\t$geneID\t$isoID\t$start\t$end\t";
    }
    if (0 == $j % 2) {   ### this is an even number
          my ($chr, $isoID, $geneID, $type, $start, $end, $strand, $cDNA_accs) = split (/\s+/, $link{$j});
          print OUTPUT "$isoID\t$start\t$end\n";
    }
}
close OUTPUT;

#chr1    retained_exon   -       6       asmbl_19        679289  679421  asmbl_20        678538  679670

### reopen the temp file
open (OUTPUT, 'temp') or die "Can not open input file: temp\n";
while (<OUTPUT>) {
               chomp;
               if ($_ =~ /AS_type/) {
                      print "$_\n";
                      next;
               }
               my @list = split (/\s+/, $_);
               if ($list[4] eq $list[7]) {  ### This means actually this is probably bug of PASA rarely see the same isoform ID should have two different isoforms ID
                      next;   
               }

               if ($_ =~ /alt_acceptor/) {
                      my @line = split (/\s+/, $_);
                      my ($A_start, $B_start, $A_end, $B_end);
                      if ( $line[2] =~ /\+/ ) {
                           $A_end = $line[6] + 1;
                           $B_end = $line[9] + 1;
                           foreach my $pos_A ( sort keys %{$exon_start{$line[4]}} ) {
                                   if ($exon_start{$line[4]}{$pos_A} =~ /$A_end/) {
                                           $A_start = $exon_end{$line[4]}{($pos_A - 1)};
                                   }
                           }
                           
                           foreach my $pos_B ( sort keys %{$exon_start{$line[7]}} ) {
                                   if ($exon_start{$line[7]}{$pos_B} =~ /$B_end/) {
                                          $B_start = $exon_end{$line[7]}{($pos_B - 1)};
                                   }
                           }
                      }
                      
                      if ( $line[2] =~ /\-/ ) {
                           $A_end = $line[5] - 1;
                           $B_end = $line[8] - 1;
                           foreach my $pos_A ( sort keys %{$exon_end{$line[4]}} ) {
                                   if ($exon_end{$line[4]}{$pos_A} =~ /$A_end/) {
                                           $A_start = $exon_start{$line[4]}{($pos_A + 1)};
                                   }
                           }

                           foreach my $pos_B ( sort keys %{$exon_end{$line[7]}} ) {
                                   if ($exon_end{$line[7]}{$pos_B} =~ /$B_end/) {
                                           $B_start = $exon_start{$line[7]}{($pos_B + 1)};
                                   }
                           }
                      }

                      my $A_dis = abs($A_end - $A_start);
                      my $B_dis = abs($B_end - $B_start);                      
                      if ($A_dis > $B_dis) {
                                   print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$A_start\t$A_end\t$line[7]\t$B_start\t$B_end\n";
                      } else {
                                   print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[7]\t$B_start\t$B_end\t$line[4]\t$A_start\t$A_end\n";
                      }
               }

               if ($_ =~ /alt_donor/) {
                      my @line = split (/\s+/, $_);
                      my ($A_start, $B_start, $A_end, $B_end);
                      if ( $line[2] =~ /\+/ ) { 
                           $A_start = $line[5] - 1;
                           $B_start = $line[8] - 1;
                           foreach my $pos_A ( sort keys %{$exon_end{$line[4]}} ) {
                                   if ($exon_end{$line[4]}{$pos_A} =~ /$A_start/) {
                                           $A_end = $exon_start{$line[4]}{($pos_A + 1)};
                                   }
                           }

                           foreach my $pos_B ( sort keys %{$exon_end{$line[7]}} ) {
                                   if ($exon_end{$line[7]}{$pos_B} =~ /$B_start/) {
                                          $B_end = $exon_start{$line[7]}{($pos_B + 1)};
                                   }
                           }
                      }

                      if ( $line[2] =~ /\-/ ) {   
                           $A_start = $line[6] + 1;
                           $B_start = $line[9] + 1;
                           foreach my $pos_A ( sort keys %{$exon_start{$line[4]}} ) {
                                   if ($exon_start{$line[4]}{$pos_A} =~ /$A_start/) {
                                           $A_end = $exon_end{$line[4]}{($pos_A - 1)};
                                   }
                           }

                           foreach my $pos_B ( sort keys %{$exon_start{$line[7]}} ) {
                                   if ($exon_start{$line[7]}{$pos_B} =~ /$B_start/) {
                                           $B_end = $exon_end{$line[7]}{($pos_B - 1)};
                                   }
                           }
                      }

                      my $A_dis = abs($A_end - $A_start);
                      my $B_dis = abs($B_end - $B_start);
                      if ($A_dis > $B_dis) {
                                   print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$A_start\t$A_end\t$line[7]\t$B_start\t$B_end\n";
                      } else {
                                   print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[7]\t$B_start\t$B_end\t$line[4]\t$A_start\t$A_end\n";
                      }
               }
                           
               if ($_ =~ /alternate_exon/) {
                      my @line = split (/\s+/, $_);
                      my ($A_start, $B_start, $A_end, $B_end);

                      if ( $line[2] =~ /\+/ ) {
                           if ($exon_start{$line[4]}{1} =~ /$line[5]/) {
                                   ### this is 5' terminal
                                   $A_start = $line[6];
                                   $B_start = $line[9];
                                   foreach my $pos_A ( sort keys %{$exon_end{$line[4]}} ) {
                                            if ($exon_end{$line[4]}{$pos_A} =~ /$A_start/) {
                                                    $A_end = $exon_start{$line[4]}{($pos_A + 1)};
                                            }
                                   }

                                   foreach my $pos_B ( sort keys %{$exon_end{$line[7]}} ) {
                                            if ($exon_end{$line[7]}{$pos_B} =~ /$B_start/) {
                                                    $B_end = $exon_start{$line[7]}{($pos_B + 1)};
                                            }
                                   }                           
                           } else { #### this is on the 3' terminal 
  
                                   $A_end = $line[5];
                                   $B_end = $line[8];
                                   foreach my $pos_A ( sort keys %{$exon_start{$line[4]}} ) {
                                            if ($exon_start{$line[4]}{$pos_A} =~ /$A_end/) {
                                                     $A_start = $exon_end{$line[4]}{($pos_A - 1)};
                                            }
                                   } 

                                   foreach my $pos_B ( sort keys %{$exon_start{$line[7]}} ) {
                                            if ($exon_start{$line[7]}{$pos_B} =~ /$B_end/) {
                                                     $B_start = $exon_end{$line[7]}{($pos_B - 1)};
                                            }
                                   }
                           }
                      }

  
                      if ( $line[2] =~ /\-/ ) { 
                           if ($exon_start{$line[4]}{1} =~ /$line[5]/) {  #### this is on the 3' terminal
                                   $A_end = $line[6];
                                   $B_end = $line[9];
                                   foreach my $pos_A ( sort keys %{$exon_end{$line[4]}} ) {
                                           if ($exon_end{$line[4]}{$pos_A} =~ /$A_end/) {
                                                     $A_start = $exon_start{$line[4]}{($pos_A + 1)};
                                           }
                                   } 
 
                                   foreach my $pos_B ( sort keys %{$exon_end{$line[7]}} ) {
                                           if ($exon_end{$line[7]}{$pos_B} =~ /$B_end/) {
                                                     $B_start = $exon_start{$line[7]}{($pos_B + 1)};
                                           }
                                   }    
                           } else { ### this is on the 5' terminal
                                   $A_start = $line[5];
                                   $B_start = $line[8];
                                   foreach my $pos_A ( sort keys %{$exon_start{$line[4]}} ) {
                                           if ($exon_start{$line[4]}{$pos_A} =~ /$A_start/) {
                                                     $A_end = $exon_end{$line[4]}{($pos_A - 1)};
                                           }
                                   }

                                   foreach my $pos_B ( sort keys %{$exon_start{$line[7]}} ) {
                                           if ($exon_start{$line[7]}{$pos_B} =~ /$B_start/) {
                                                     $B_end = $exon_end{$line[7]}{($pos_B - 1)};
                                           }
                                   }
                           }
                      }
                      
                      my $A_dis = abs($A_end - $A_start);
                      my $B_dis = abs($B_end - $B_start);
                      if ($A_dis > $B_dis) {
                                   print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$A_start\t$A_end\t$line[7]\t$B_start\t$B_end\n";
                      } else {
                                   print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[7]\t$B_start\t$B_end\t$line[4]\t$A_start\t$A_end\n";
                      }
               }
               
 
               if ($_ =~ /retained_intron/) {
                      my @line = split (/\s+/, $_);
                      ### convert intron based location to exon based location
                      my $skip_left  = $line[8] - 1;
                      my $skip_right = $line[9] + 1;
                      if ( $line[2] =~ /\+/ ) { 
                           print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\tNA\tNA\t$line[7]\t$skip_left\t$skip_right\n";
                      }
                      if ( $line[2] =~ /\-/ ) {
                           print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\tNA\tNA\t$line[7]\t$skip_right\t$skip_left\n";
                      }
               }

               if ($_ =~ /retained_exon/) {  #### change to skipped exon, the first isofrom is the skip isoform
                      my @line = split (/\s+/, $_); 
                      ### convert intron based location to exon based location
                      my $skip_left  = $line[8] - 1;
                      my $skip_right = $line[9] + 1;
                      if ( $line[2] =~ /\+/ ) { 
                           print "$line[0]\tskipped_exon\t$line[2]\t$line[3]\t$line[7]\t$skip_left\t$skip_right\t$line[4]\t$line[5]\t$line[6]\n";
                      }
                      if ( $line[2] =~ /\-/ ) { 
                           print "$line[0]\tskipped_exon\t$line[2]\t$line[3]\t$line[7]\t$skip_right\t$skip_left\t$line[4]\t$line[6]\t$line[5]\n";
                      }
               }
}
close OUTPUT;
        
unlink "temp";


















