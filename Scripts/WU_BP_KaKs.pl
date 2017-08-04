#!/sw/bin/perl5.12.3 -w
use strict;
use BPlitenew;

### output of this scripts will produce a axt format file to be used in the KaKs_Calculator2.0

my $usage = '
WU_BP_KaKs.pl <BLAST REPORT> <conserved AS file> <fasta file for junction seq>
';
die $usage unless (-e $ARGV[2]);

# conserved_IntronR_synteny_5Grass.csv 
# IntronR_1983_2971	3	NA	NA	LOC_Os09g27010:IntronR_Chr9_-_NA_NA_16424672_16424550;	Bradi4g31430:IntronR_Bd4_-_NA_NA_37274066_37273979;	Si029783m.g:IntronR_scaffold_2_-_NA_NA_32210082_32209994;

my %conserved_AS;
open (AS, "$ARGV[1]") or die "Can not open $ARGV[1]\n";
while (<AS>) {
	  chomp;
          next if $_ =~ /^ClusterID/;
          my ($ClusterID, $Num_Of_Species, $Maize, $Sorghum, $Rice, $Brachypodium, $Millet) = split (/\s+/, $_);
	  my @array_maize = split (/\;/, $Maize);
	  foreach my $i (@array_maize) {
		next if $Maize =~ /^NA/;
		if ($i !~ /"NA"/) {
			$conserved_AS{$i} = $ClusterID;
		}
	  }
	  my @array_sorghum = split (/\;/, $Sorghum);
          foreach my $i (@array_sorghum) {
                next if $Sorghum =~ /^NA/;
                if ($i !~ /"NA"/) {
                        $conserved_AS{$i} = $ClusterID;
                }  
	  }
	  my @array_rice = split (/\;/, $Rice);
          foreach my $i (@array_rice) {
                next if $Rice =~ /^NA/;
                if ($i !~ /"NA"/) {
                        $conserved_AS{$i} = $ClusterID;
                }
          }
	  my @array_brachy = split (/\;/, $Brachypodium);
          foreach my $i (@array_brachy) {
                next if $Brachypodium =~ /^NA/;
                if ($i !~ /"NA"/) {
                        $conserved_AS{$i} = $ClusterID;
                }
          }
	  my @array_millet = split (/\;/, $Millet);
          foreach my $i (@array_millet) {
                next if $Millet =~ /^NA/;
                if ($i !~ /"NA"/) {
                        $conserved_AS{$i} = $ClusterID;
                }
          }
}
close AS;

my %junc_seq;
my ($seq, $name);
open (FASTA, "$ARGV[2]") or die "Can not open $ARGV[2]\n";
while (<FASTA>) {
          chomp;
	  if (/^>(\S+)/) {
                if (defined $seq) {
			#print "$name\n"; 
                        $junc_seq{$name} = $seq;
                }       
                $name = $1;
                $seq = "";
        } else {
                $seq .= $_; 
        }
}
$junc_seq{$name} = $seq;
close FASTA;


open (BLAST, "$ARGV[0]") or die "Can not open $ARGV[0]\n";
#for blastx, query alignment coordinates in NT, sub coords in aa,  Alignment length based on AA identities based on AAs
#take only top hit, hash HSPs, and process to get longest alignment
my $multiple_report = new BPlite::Multi(\*BLAST );
 #print "mult rep = $multiple_report\n";
my %output; 
 #for each BLAST report.....
 QUERY:	while(my $blast = $multiple_report->nextReport) {
	 my $query = $blast->query;
	 #print "#$query\n";
	 $query =~ s/se q/seq/;
	 $query =~ s/s eq/seq/;
	 $query =~ s/seq 1/seq1/;
	 $query =~ s/seq 2/seq2/;
	 $query =~ s/ seq1/seq1/;
	 $query =~ s/ seq2/seq2/;
	 #print "#$query\n";
	 my ($qname, $qlength) = $query =~ /(\S+).+\((\d+,*\d+)\sletters/;
	 #print "$qname\t$qlength\n";
	 my @q_array = split (/\:/, $qname);
	 my $qname_partial = $q_array[0].":".$q_array[1];
	 next QUERY if (not defined $conserved_AS{$qname_partial}); 
	 SUBJECT: while(my $subject = $blast->nextSbjct) {
         #print "query name = $qname\n";
		 #next QUERY if $subject =~ /transpos|transcriptase|retro|gag|pol|gag-pol|polyprotein|mudr|rire|copia|poson|olyprotein|MITE|IS-element/i;
		 #TIGR ATH specific
		 my ($sname) = $subject->name =~ />(\S+)/;
		 my @s_array = split (/\:/, $sname);
		 my $sname_partial = $s_array[0].":".$s_array[1];
		 next SUBJECT if ( $q_array[0] =~ /^Si/ && $s_array[0] =~ /^Si/);
		 next SUBJECT if ( $q_array[0] =~ /^LOC/ && $s_array[0] =~ /^LOC/);
		 next SUBJECT if ( $q_array[0] =~ /^Sobic/ && $s_array[0] =~ /^Sobic/);
		 next SUBJECT if ( $q_array[0] =~ /^GRMZM/ && $s_array[0] =~ /^GRMZM/);
		 next SUBJECT if ( $q_array[0] =~ /^Bradi/ && $s_array[0] =~ /^Bradi/);
		 next SUBJECT if ($q_array[0] =~ $s_array[0]);
		 next SUBJECT if (not defined $conserved_AS{$sname_partial});
		 my $slength = $subject->length;
		 while(my $hsp = $subject->nextHSP) {
		 	my $score =$hsp->score;
			my $bits = $hsp->bits;
			my $id = $hsp->percent;
			my $pval = $hsp->P;
			my $match = $hsp->match;
			my $pos = $hsp->positive;
			my $alig_len = $hsp->length;
			my $qb = $hsp->queryBegin;	  
			my $qe = $hsp->queryEnd; 	  
			my $sb = $hsp->sbjctBegin;	  
			my $se = $hsp->sbjctEnd; 	  
			my $q_alig = $hsp->queryAlignment; 
			my $sub_alig = $hsp->sbjctAlignment; 
			my $alig_string = $hsp->alignmentString;
			my $q_gaps = $hsp->queryGaps;	  
			my $s_gaps = $hsp->sbjctGaps;	  
			my $strand;
			my $alig_p;
			if ($qb <=$qe){
				$alig_p = $qe - $qb;
				$strand = "+";
			} elsif ($qe <=$qb){
				 $alig_p = $qb - $qe;
				 $strand = "-";
			}
			my $qdis = abs($qe - $qb);
			my $sdis = abs($se - $sb);
			#my $allowed_overhang = (0.10 * $alig_len);
			#my $cons = ($pos/$alig_len)*100;
		 	next if ($qname_partial eq $sname_partial);
			#$alig_p =~ s/\,//g; $qlength =~ s/\,//g;
			#my $alig_test = $alig_p/$qlength;
			# ($pval <= 1e-100) 
			#next unless ($id >= 95);# && (abs($se-$sb)/$slength >= 0.80);# && ($alig_len > 450);# || ($alig_test >= 0.90)) && ($alig_len > 450); 
			#my $xyFLAG = XY($qlength, $slength, $qb, $qe, $sb, $se, $strand, $allowed_overhang); 
			#next unless ($qlength >= 1000) && ($slength>=1000);# && ($xyFLAG == 0); #&& ($qlength <= 5000) && ($slength<=5000);
		        if ($qdis >= 120 && $sdis >= 120 && not defined $output{$qname}{$sname}) {
				#print "#####$qname\t*******$sname\t$alig_len\n";
				my $nt_subject;
				my $nt_query;
				my $qpos_nt = $qb - 1;
				my $spos_nt = $sb - 1;
				for (my $i=1; $i<=$alig_len; $i++) {
					my $q_AA = substr($q_alig, $i-1, 1);
					if ($q_AA !~ /\*/ && $q_AA !~ /\-/) {
						$nt_query = $nt_query.substr($junc_seq{$qname}, $qpos_nt, 3);
						#print "**************************$nt_subject\t#####$qname\t*******$sname\t$alig_len\n";
						$qpos_nt = $qpos_nt + 3;
					}
					if ($q_AA =~ /\*/) { ### stop codon
						$nt_query = $nt_query."---";
						$qpos_nt = $qpos_nt + 3;
					}
					if ($q_AA =~ /\-/) { ### match blank space
                                                $nt_query = $nt_query."---";
                                                $qpos_nt = $qpos_nt + 0;
                                        }
				}
				for (my $i=1; $i<=$alig_len; $i++) {
                                        my $s_AA = substr($sub_alig, $i-1, 1);
                                        if ($s_AA !~ /\*/ && $s_AA !~ /\-/) {
						$nt_subject = $nt_subject.substr($junc_seq{$sname}, $spos_nt, 3);
						$spos_nt = $spos_nt + 3;
					}
					if ($s_AA =~ /\*/) {
						$nt_subject = $nt_subject."---";
                                                $spos_nt = $spos_nt + 3;
                                        }
                                        if ($s_AA =~ /\-/) {
						$nt_subject = $nt_subject."---";
                                                $spos_nt = $spos_nt + 0;
                                        }
                                }
				print "$qname"."@"."$sname\n";
				print "$nt_query\n";
				print "$nt_subject\n";
				print "\n";
				$output{$qname}{$sname}++;
				$output{$sname}{$qname}++;
			}
		#	print "$qname\t$qlength\t$sname\t$slength\t$score\t$pval\t$id\t$cons\t$alig_len\t$qb\t$qe\t$sb\t$se\t$strand\n" ;#if ( $id >= 65 && $alig_len >= 150);
	#next QUERY;
		 }
		 
		 
     }
    
 }

#QUERYNAME, QUERYLEN, SUBNAME, SUBLEN, SCORE, PVAL, ID, CONS, ALIG_LEN, QB, QE, SB, SE

		
__END__

###################
#column definitions
###################

#################################
# parseblast column definitions #
#################################

 $hsp->score;
 $hsp->bits;
 $hsp->percent;
 $hsp->P;
 $hsp->match;
 $hsp->positive;
 $hsp->length;
 $hsp->queryBegin;      $hsp->qb;
 $hsp->queryEnd;        $hsp->qe;
 $hsp->sbjctBegin;      $hsp->sb;
 $hsp->sbjctEnd;        $hsp->se;
 $hsp->queryAlignment;  $hsp->qa;
 $hsp->sbjctAlignment;  $hsp->sa;
 $hsp->alignmentString; $hsp->as;
 $hsp->queryGaps;       $hsp->qg;
 $hsp->sbjctGaps;       $hsp->sg;









 1: query begin
 2: query end
 3: query name
 4: sbjct begin
 5: sbjct end
 6: sbjct name
 7: raw score
 8: bits (normalized score)
 9: E-value
10: P-value (really, also the E value)
11: percent identity
12: number of matches
13: number of positive scores (similarities)
14: length of alignment
15: length of query
16: length of sbjct
17: number of gaps in query alignment
18: number of gaps in sbjct alignment

END

