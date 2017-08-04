#!/usr/bin/perl -w
use strict;

### Modify Oct23, 2015, this version I separate two sides of junction instead of connect them into one piece of sequence
### 1) This scripts take AS event file, GTF file, genome fasta file and cut-off length to identify 
### 2) the junction sequence on the two sides of exons, the length is up to 300bp or just the length of the exon
### 3) if there are multiple isoforms support the same event, pick up the longest franking sequence

my $usage = "$0 <Final_PASA_OUTPUT_pasa_wmei_B73_March18th_2015.indiv_splice_labels_and_coords.dat> <Final_FilterSingleTranscriptLoci_B73_entropyscore2_IRfilter_fpkm1frac0.05.gtf> <genome fasta sequence> <length cut off>";
my ($PASA_AS, $PASA_GTF, $genome_fa, $Length_cutoff) = @ARGV;
my (%chr_seq, %exon_start, %exon_end, %list_transID);

#Chr     AS_type strand  Gene_ID Isoform_A       A_start A_end   Isoform_B       B_start B_end
#chr1    skipped_exon    -       GRMZM2G059745   asmbl_27        93386   92622   asmbl_35        93242   92905
#chr1    alt_acceptor    -       GRMZM2G059745   asmbl_27        93386   92622   asmbl_36        93386   93242
#chr1    retained_intron -       GRMZM2G059745   asmbl_36        NA      NA      asmbl_35        92905   92622
#chr1    retained_intron -       GRMZM5G833153   asmbl_41        NA      NA      asmbl_43        144990  143713
#chr1    retained_intron -       GRMZM5G833153   asmbl_43        NA      NA      asmbl_41        145219  145069
#chr1    alt_acceptor    -       GRMZM5G833153   asmbl_41        145219  145069  asmbl_47        145219  145161
#chr1    alt_acceptor    -       GRMZM5G833153   asmbl_41        145219  145069  asmbl_48        145219  145178

#chr1    assembler       transcript      84675   93563   .       -       .       gene_id "GRMZM2G059745"; transcript_id "align_id:1142745|asmbl_27";
#chr1    assembler       exon    84675   85000   .       -       .       gene_id "GRMZM2G059745"; transcript_id "align_id:1142745|asmbl_27";
#chr1    assembler       exon    86014   86169   .       -       .       gene_id "GRMZM2G059745"; transcript_id "align_id:1142745|asmbl_27";
#chr1    assembler       exon    90685   90798   .       -       .       gene_id "GRMZM2G059745"; transcript_id "align_id:1142745|asmbl_27";
#chr1    assembler       exon    92521   92622   .       -       .       gene_id "GRMZM2G059745"; transcript_id "align_id:1142745|asmbl_27";
#chr1    assembler       exon    93386   93563   .       -       .       gene_id "GRMZM2G059745"; transcript_id "align_id:1142745|asmbl_27";

my ($transID, $count);

### Read the PASA GTF exon table
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
               $count++;
          }
}     
close GTF;

my $IR_output = "IntronR_seq1seq2_$Length_cutoff"."_".$PASA_AS.".fasta";
my $Alt_acceptor_output = "AltA_seq1seq2_$Length_cutoff"."_".$PASA_AS.".fasta";
my $Alt_donor_output = "AltD_seq1seq2_$Length_cutoff"."_".$PASA_AS.".fasta";
my $Exon_skip_output = "ExonS_seq1seq2_$Length_cutoff"."_".$PASA_AS.".fasta";
my $ATE_output = "ATE_seq1seq2_$Length_cutoff"."_".$PASA_AS.".fasta";
my $output_file = "AddFrankExonSeq_seq1seq2_$Length_cutoff"."_".$PASA_AS;

open (OUTPUT, ">$output_file"); 
open (IntronR, ">$IR_output");
open (AltA, ">$Alt_acceptor_output");
open (AltD, ">$Alt_donor_output");
open (ExonS, ">$Exon_skip_output");
open (ATE, ">$ATE_output");

my ($seq, $chr);

### read the genome file
open (IN, $genome_fa) or die "Can not open input file: $genome_fa\n";
while (<IN>) {
        chomp;
        if (/^>(\S+)/) {
                if (defined $seq) {
                        $chr_seq{$chr} = $seq;
                }
                $chr = $1;
                $seq = "";
        } else {
                $seq .= $_;
        }
}
$chr_seq{$chr} = $seq;
close IN;


open (AS, $PASA_AS) or die "Can not open input file: $PASA_AS\n";
while (<AS>) {
               chomp;
               if ($_ =~ /AS_type/) {
                      print OUTPUT "$_\tflank_exon_seq\n";
                      next;
               }
               my ($upstream_frank_seq, $downstream_frank_seq);
               my ($Chr, $AS_type, $strand, $Gene_ID, $Isoform_A, $A_start, $A_end, $Isoform_B, $B_start, $B_end) = split (/\s+/, $_);
               my ($Isoform_ID, $Junc_start, $Junc_end);

               if ($AS_type =~ /[skipped_exon|alt_acceptor|alt_donor|alternate_exon]/) {
		      $Isoform_ID = $Isoform_A;
                      $Junc_start = $A_start;
                      $Junc_end   = $A_end;
	       }
            
               if ($AS_type =~ /retained_intron/) {
                      $Isoform_ID = $Isoform_B;
                      $Junc_start = $B_start;
                      $Junc_end   = $B_end;
               }
                      
               if ($strand =~ /\+/) {
                   	      foreach my $element (sort keys %{$exon_start{$Isoform_ID}}) {
					if ($exon_start{$Isoform_ID}{$element} == $Junc_end) {
                                                if ( ($exon_end{$Isoform_ID}{$element-1} - $exon_start{$Isoform_ID}{$element-1} + 1) <= $Length_cutoff) {
                                                       my $start = $exon_start{$Isoform_ID}{$element-1} - 1;
						       my $distance = $exon_end{$Isoform_ID}{$element-1} - $exon_start{$Isoform_ID}{$element-1} + 1;
                                                       $upstream_frank_seq = substr($chr_seq{$Chr}, $start, $distance);
                                                }

                                                if ( ($exon_end{$Isoform_ID}{$element-1} - $exon_start{$Isoform_ID}{$element-1} + 1) > $Length_cutoff) {
                                                       my $start = $exon_start{$Isoform_ID}{$element-1} - 1;
                                                       $upstream_frank_seq = substr($chr_seq{$Chr}, $start, $Length_cutoff);
                                                }
     
                                                if ( ($exon_end{$Isoform_ID}{$element} - $Junc_end + 1) <= $Length_cutoff) {
                                                       my $start = $Junc_end - 1;
                                                       my $distance = $exon_end{$Isoform_ID}{$element} - $Junc_end + 1;
                                                       $downstream_frank_seq = substr($chr_seq{$Chr}, $start, $distance);
                                                }

						if ( ($exon_end{$Isoform_ID}{$element} - $Junc_end + 1) > $Length_cutoff) {
                                                       my $start = $Junc_end - 1;
                                                       $downstream_frank_seq = substr($chr_seq{$Chr}, $start, $Length_cutoff);
                                                }
 
                                                my $combined_seq = $upstream_frank_seq.$downstream_frank_seq;
                                                print OUTPUT "$_\t$upstream_frank_seq#$downstream_frank_seq#$combined_seq\n";        
						#print "$exon_start{$Isoform_ID}{$element-1}\t$exon_end{$Isoform_ID}{$element-1}\t$Junc_end\t$exon_end{$Isoform_ID}{$element}\t$_\t$combined_seq\n";
                                        }
                              }
               }
  
               if ($strand =~ /\-/) {
                              foreach my $element (sort keys %{$exon_start{$Isoform_ID}}) {
                                        if ($exon_start{$Isoform_ID}{$element} == $Junc_start) {
						my ($revcom_upstream_frank_seq, $revcom_downstream_frank_seq);
                                                if ( ($exon_end{$Isoform_ID}{$element-1} - $exon_start{$Isoform_ID}{$element-1} + 1) <= $Length_cutoff) {
                                                       my $start = $exon_start{$Isoform_ID}{$element-1} - 1;
                                                       my $distance = $exon_end{$Isoform_ID}{$element-1} - $exon_start{$Isoform_ID}{$element-1} + 1;
                                                       $upstream_frank_seq = substr($chr_seq{$Chr}, $start, $distance);
						       ### revcom the seq
						       $revcom_upstream_frank_seq = reverse $upstream_frank_seq;
						       $revcom_upstream_frank_seq =~ tr/ACGTacgt/TGCAtgca/;
                                                }

                                                if ( ($exon_end{$Isoform_ID}{$element-1} - $exon_start{$Isoform_ID}{$element-1} + 1) > $Length_cutoff) {
                                                       my $start = $exon_start{$Isoform_ID}{$element-1} - 1;
                                                       $upstream_frank_seq = substr($chr_seq{$Chr}, $start, $Length_cutoff);
						       $revcom_upstream_frank_seq = reverse $upstream_frank_seq;
                                                       $revcom_upstream_frank_seq =~ tr/ACGTacgt/TGCAtgca/;
                                                }

                                                if ( ($exon_end{$Isoform_ID}{$element} - $Junc_start + 1) <= $Length_cutoff) {
                                                       my $start = $Junc_start - 1;
                                                       my $distance = $exon_end{$Isoform_ID}{$element} - $Junc_start + 1;
                                                       $downstream_frank_seq = substr($chr_seq{$Chr}, $start, $distance);
						       $revcom_downstream_frank_seq = reverse $downstream_frank_seq;
                                                       $revcom_downstream_frank_seq =~ tr/ACGTacgt/TGCAtgca/;	
                                                }

                                                if ( ($exon_end{$Isoform_ID}{$element} - $Junc_start + 1) > $Length_cutoff) {
                                                       my $start = $Junc_start - 1;
                                                       $downstream_frank_seq = substr($chr_seq{$Chr}, $start, $Length_cutoff);
						       $revcom_downstream_frank_seq = reverse $downstream_frank_seq;
                                                       $revcom_downstream_frank_seq =~ tr/ACGTacgt/TGCAtgca/;	
                                                }

                                                my $combined_seq = $revcom_downstream_frank_seq.$revcom_upstream_frank_seq;
                                                print OUTPUT "$_\t$revcom_downstream_frank_seq#$revcom_upstream_frank_seq#$combined_seq\n";
                                        }
                               }
              }
}
close AS;
close OUTPUT;

my (%alt_acceptor_list, %alt_donor_list, %alternate_exon_list, %skipped_exon_list, %retained_intron_list, %gene_link);

open (INPUT, $output_file) or die "Can not open input file: $output_file\n";
while (<INPUT>) {
		chomp;
		next if ($_ =~ /AS_type/);
                my ($Chr, $AS_type, $strand, $Gene_ID, $Isoform_A, $A_start, $A_end, $Isoform_B, $B_start, $B_end, $flank_exon_seq) = split (/\s+/, $_);
	        ### because the way I wrote the $flank_exon_seq = $seq1#$seq2#$combined_seq, the longer of the $flank_exon_seq is equal to longer of $combined_seq;
		my $seq_length = length($flank_exon_seq);
		#print "$_\t$flank_exon_seq\n";
                my $AS_event_ID;
                if ($AS_type =~ /alt_acceptor/) {
			$AS_event_ID = "AltA"."_".$Chr."_".$strand."_".$A_start."_".$A_end."_".$B_start."_".$B_end.":";
	                if (not defined $gene_link{$AS_event_ID}) {
            			$gene_link{$AS_event_ID} = $Gene_ID;
			}
			if (not defined $alt_acceptor_list{$AS_event_ID}) {
                                $alt_acceptor_list{$AS_event_ID} = $flank_exon_seq;
                        }

                        if ( defined $alt_acceptor_list{$AS_event_ID} ) {
 				my $old_seq_length = length($alt_acceptor_list{$AS_event_ID});
                                if ($old_seq_length < $seq_length) {
					$alt_acceptor_list{$AS_event_ID} = $flank_exon_seq;
				}
			}
		}

                if ($AS_type =~ /alt_donor/) {
                        $AS_event_ID = "AltD"."_".$Chr."_".$strand."_".$A_start."_".$A_end."_".$B_start."_".$B_end.":";
                        if (not defined $gene_link{$AS_event_ID}) {
				$gene_link{$AS_event_ID} = $Gene_ID;
			}
			if (not defined $alt_donor_list{$AS_event_ID}) {
                                $alt_donor_list{$AS_event_ID} = $flank_exon_seq;
                        }
             
                        if ( defined $alt_donor_list{$AS_event_ID} ) {
                                my $old_seq_length = length($alt_donor_list{$AS_event_ID});
                                if ($old_seq_length < $seq_length) {
                                        $alt_donor_list{$AS_event_ID} = $flank_exon_seq;
                                }
                        }
                }

		if ($AS_type =~ /alternate_exon/) {
                        $AS_event_ID = "ATE"."_".$Chr."_".$strand."_".$A_start."_".$A_end."_".$B_start."_".$B_end.":";
                        if (not defined $gene_link{$AS_event_ID}) {
				$gene_link{$AS_event_ID} = $Gene_ID;
			}
			if (not defined $alternate_exon_list{$AS_event_ID}) {
                                $alternate_exon_list{$AS_event_ID} = $flank_exon_seq;
                        }
                        if ( defined $alternate_exon_list{$AS_event_ID} ) {
                                my $old_seq_length = length($alternate_exon_list{$AS_event_ID});
                                if ($old_seq_length < $seq_length) {
                                        $alternate_exon_list{$AS_event_ID} = $flank_exon_seq;
                                }
                        }
                }

		if ($AS_type =~ /skipped_exon/) {
                        $AS_event_ID = "ExonS"."_".$Chr."_".$strand."_".$A_start."_".$A_end."_NA_NA:";
                        if (not defined $gene_link{$AS_event_ID}) {
				$gene_link{$AS_event_ID} = $Gene_ID;
			}
			if (not defined $skipped_exon_list{$AS_event_ID}) {
                                $skipped_exon_list{$AS_event_ID} = $flank_exon_seq;
                        }
                        
                        if ( defined $skipped_exon_list{$AS_event_ID} ) {
                                my $old_seq_length = length($skipped_exon_list{$AS_event_ID});
                                if ($old_seq_length < $seq_length) {
                                        $skipped_exon_list{$AS_event_ID} = $flank_exon_seq;
                                }
                        }
                }

		if ($AS_type =~ /retained_intron/) {
                        $AS_event_ID = "IntronR"."_".$Chr."_".$strand."_NA_NA_".$B_start."_".$B_end.":";
                        if (not defined $gene_link{$AS_event_ID}) {
				$gene_link{$AS_event_ID} = $Gene_ID;
			}
			if (not defined $retained_intron_list{$AS_event_ID}) {
                                $retained_intron_list{$AS_event_ID} = $flank_exon_seq;
                        }
                        
                        if ( defined $retained_intron_list{$AS_event_ID} ) {
                                my $old_seq_length = length($retained_intron_list{$AS_event_ID});
                                if ($old_seq_length < $seq_length) {
                                        $retained_intron_list{$AS_event_ID} = $flank_exon_seq;
                                }
                        }
                }		
}
close INPUT;

foreach my $key1 (keys %alt_acceptor_list) {
	my $geneid = $gene_link{$key1};
	my ($seq1, $seq2, $combined_seq) = split (/\#/, $alt_acceptor_list{$key1});
	my $len_seq1 = length($seq1);
	my $len_seq2 = length($seq2);
	print AltA ">$geneid:$key1$len_seq1:seq1\n$seq1\n";
	print AltA ">$geneid:$key1$len_seq2:seq2\n$seq2\n";
}

foreach my $key2 (keys %alt_donor_list) {
	my $geneid = $gene_link{$key2};
	my ($seq1, $seq2, $combined_seq) = split (/\#/, $alt_donor_list{$key2});
        my $len_seq1 = length($seq1);
        my $len_seq2 = length($seq2);
        print AltD ">$geneid:$key2$len_seq1:seq1\n$seq1\n";
        print AltD ">$geneid:$key2$len_seq2:seq2\n$seq2\n";
}

foreach my $key3 (keys %alternate_exon_list) {
	my $geneid = $gene_link{$key3};
	my ($seq1, $seq2, $combined_seq) = split (/\#/, $alternate_exon_list{$key3});
        my $len_seq1 = length($seq1);
        my $len_seq2 = length($seq2);
        print ATE ">$geneid:$key3$len_seq1:seq1\n$seq1\n";
        print ATE ">$geneid:$key3$len_seq2:seq2\n$seq2\n";
}

foreach my $key4 (keys %skipped_exon_list) {
	my $geneid = $gene_link{$key4};
	my ($seq1, $seq2, $combined_seq) = split (/\#/, $skipped_exon_list{$key4});
        my $len_seq1 = length($seq1);
        my $len_seq2 = length($seq2);
        print ExonS ">$geneid:$key4$len_seq1:seq1\n$seq1\n";
        print ExonS ">$geneid:$key4$len_seq2:seq2\n$seq2\n";
}

foreach my $key5 (keys %retained_intron_list) {
	my $geneid = $gene_link{$key5};
	my ($seq1, $seq2, $combined_seq) = split (/\#/, $retained_intron_list{$key5});
        my $len_seq1 = length($seq1);
        my $len_seq2 = length($seq2);
        print IntronR ">$geneid:$key5$len_seq1:seq1\n$seq1\n";
        print IntronR ">$geneid:$key5$len_seq2:seq2\n$seq2\n";
}


