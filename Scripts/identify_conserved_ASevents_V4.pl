#!/usr/bin/perl -w
use strict;
use lib "/ufrc/barbazuk/wmei/PROJECT/FROM_SCRATCH/src/Module";
use Utils qw(:all);

### this scripts takes the tblastx output file and choose the output E-value <= 1e-5
### require at least both two side has 30bp
### generate the conserved clustering table for each AS event type
### NO E-value cutoff
### remove W22 genome from this analysis (Sep6, 2016)

my $usage = "$0 <tblastx output> <OrthologousGroups.txt> <AS_event_type>";
die $usage unless $#ARGV == 2;
my ($tblastx_output, $OrthoGroup, $AS_event_type) = @ARGV;
my (%clusterID, %keep);

# Fields: qid   sid     E       N       Sprime  S       alignlen        nident  npos    nmism   pcident pcpos   qgaps   qgaplen sgaps   sgaplen qframe  qstart  qend    sframe  sstart  send
# # EXIT: [Zm00004a002763:ExonS_Chr10_-_130300018_130292269_NA_NA:30:seq2]: 0
# # Query: Zm00004a022640:ExonS_Chr1_+_154913818_154914812_NA_NA:108:seq1
# Zm00004a022640:ExonS_Chr1_+_154913818_154914812_NA_NA:108:seq1  Zm00004a022640:ExonS_Chr1_+_154913818_154914812_NA_NA:108:seq1  3.6e-17 1       80.39   214     35      35      35      0       100.00  100.00  0       0       0       0
#        +3      3       107     +3      3       107

#GRMZM|AC|AY|EF|AF

my $i = 0;
my (%array_num, %match, %match_list, %ortho_link);
#step1: read the blast output and generate the cluster for the query and subject
open (ORTHO, $OrthoGroup) or die "Can not open input file: $OrthoGroup\n";
while (<ORTHO>) {
	  chomp;
	  my @ortholist = split(/\s+/, $_);
	  my $orthosize = @ortholist;
	  #print "$ortholist[0]\t$ortholist[$orthosize-1]\n";
	  foreach my $gene (@ortholist) {
		if ($gene =~ /B73/) {
			if ($gene =~ /B73_(GRMZM\S+)_T\d+/) {
				$gene = $1;
				$ortho_link{$gene} = $ortholist[0];
			}
			if ($gene =~ /B73_([AC|AY|EF|AY]\S+)\|PACid/) {
				$gene = $1;
				$gene =~ s/FGT/FG/;
				$ortho_link{$gene} = $ortholist[0];
			}
		}
		#if ($gene =~ /W22/) {
		#	if ($gene =~ /W22_(Zm\S+)_T/) {
		#		$gene = $1;
		#		$ortho_link{$gene} = $ortholist[0];
		#	}
		#}
		if ($gene =~ /SB/) {
			if ($gene =~ /SB_(Sobic\S+)\.\d+\.p/) {
				$gene = $1;
				$ortho_link{$gene} = $ortholist[0];
			}
		}
		if ($gene =~ /OSA/) {
			if ($gene =~ /OSA_(LOC_\S+)\.\d+/) {
				$gene = $1;
				$ortho_link{$gene} = $ortholist[0];
			}
		}
		if ($gene =~ /BD/) {
			if ($gene =~ /BD_(Bradi\S+)\.\d+\.p/) {
				$gene = $1;
				$ortho_link{$gene} = $ortholist[0];
			}
		}
		if ($gene =~ /SIT/) {
			if ($gene =~ /SIT_(Si\S+)/) {
				$gene = $1.".g";
				$ortho_link{$gene} = $ortholist[0];
			}
		}
		if ($gene =~ /MUSA/) {
			if ($gene =~ /MUSA_(GSMUA_\S+_\S+)/) {
				$gene = $1;
				$gene =~ s/P/G/;
				$ortho_link{$gene} = $ortholist[0];
			}
		}
		if ($gene =~ /EG/) {
			if ($gene =~ /EG_(p5_\S+)/) {
				$gene = $1;
				$ortho_link{$gene} = $ortholist[0];
			}
		}
		if ($gene =~ /ATH/) {
			if ($gene =~ /ATH_(AT\S+)\.\d+/) {
				$gene = $1;
				$ortho_link{$gene} = $ortholist[0];
			}
		}
		if ($gene =~ /AMBO/) {
			if ($gene =~ /AMBO_evm_27\.model\.(AmTr\S+)/) {
				$gene = $1;
				$ortho_link{$gene} = $ortholist[0];
			}
		}        	
	  }
}
close ORTHO;

open (BLAST, $tblastx_output) or die "Can not open input file: $tblastx_output\n";
while (<BLAST>) {
			chomp;
			next if /^#/; ### remove the header line
			my ($qid, $sid, $E, $N, $Sprime, $S, $alignlen, $nident, $npos, $nmism, $pcident, $pcpos, $qgaps, $qgaplen, $sgaps, $sgaplen, $qframe, $qstart, $qend, $sframe, $sstart, $send) = split (/\s+/, $_);
			### remove the blast hit with the same gene ID
			my ($genename_qid, $location_qid, $length_qid, $qseq_type) = split(/\:/, $qid);
			my ($genename_sid, $location_sid, $length_sid, $sseq_type) = split(/\:/, $sid);
			### the junction length at least 30bp for both query and subject
	  		if ($length_qid < 30 || $length_sid < 30) {
	  			next;
	  		}
			### skip query or subject are W22
			if ( $genename_qid =~ /Zm/ || $genename_sid =~ /Zm/ ) {
				next;
			}
	  		### if both query and subject more than 100bp, then at least E-value 1e-5
	  		### and if any query and subject is arabidopsis or amborella, then does not require E-value considering the long phylogenetic distance
	  		if ( $length_qid > 100 && $length_sid > 100 && $E > 1e-5 && $genename_qid !~ /[AT|AmTr]/ && $genename_sid !~ /[AT|AmTr]/) {
	  			next;
	  		}
	  		### query and subject can not be the same gene
	  		if ($genename_qid eq $genename_sid) {
	  			next;
	  		}
	  		### the gene name not in the cluster
	  		if ( (not defined $ortho_link{$genename_qid}) || (not defined $ortho_link{$genename_sid}) ) {
				next;
			}
			### query and subject belong to different cluster ID
			if ($ortho_link{$genename_qid} ne $ortho_link{$genename_sid}) {
				#print "I found it\t$genename_qid\t$ortho_link{$genename_qid}\t$genename_sid\t$ortho_link{$genename_sid}\n";
				next;
			}
			### if this match pair has been identified earlier, then move onto the next
			if ( defined $keep{$genename_qid}{$location_qid}{$genename_sid}{$location_sid} ) {
				next;
			}
			if (not defined $array_num{$genename_qid.$location_qid}) {
				$i++;
				$array_num{$genename_qid.$location_qid}++;
				$match_list{$i} = [];
				$clusterID{$i} = $ortho_link{$genename_qid};
				%keep = ();
	  		}
			#print "$location_qid\t$location_sid\t$qseq_type\t$sseq_type\n";
			$match{$location_qid}{$location_sid}{$qseq_type}{$sseq_type}++;
			if ( ( (defined $match{$location_qid}{$location_sid}{seq1}{seq1}) && (defined $match{$location_qid}{$location_sid}{seq2}{seq2}) ) || ( (defined $match{$location_qid}{$location_sid}{seq1}{seq2}) && (defined $match{$location_qid}{$location_sid}{seq2}{seq1}) ) ) {
			### only push when the value is not exist in the hash
				$qid = $genename_qid.":".$location_qid;
				$sid = $genename_sid.":".$location_sid;
				push (@{ $match_list{$i} }, $qid) unless grep{$_ eq $qid} @{ $match_list{$i} };
				push (@{ $match_list{$i} }, $sid) unless grep{$_ eq $sid} @{ $match_list{$i} };
				$keep{$genename_qid}{$location_qid}{$genename_sid}{$location_sid}++;
				$keep{$genename_sid}{$location_sid}{$genename_qid}{$location_qid}++;
			}
}     
close BLAST;

#step2: go over each list and merge overlap clustering

for (my $j = 1; $j < $i; $j++) {
	for (my $g = $j + 1; $g <= $i; $g++) {
		#print "$j\t$g\n";
		my $arrSize = scalar intersect(@{ $match_list{$j} }, @{ $match_list{$g} });
		#print "$arrSize\t$j\t$g\n";	
		if ( $arrSize > 0 ) {
			#print "I fount it\t$j!!!!!!!$g\n";
			@{ $match_list{$j} } = unique( @{ $match_list{$j} }, @{ $match_list{$g} });
			undef @{ $match_list{$g} };	
			#undef $match_list{$g};
		}
	}
}


### do the clustering twice
for (my $k = 1; $k < $i; $k++) {
        for (my $l = $k + 1; $l <= $i; $l++) {
                my $arrSize = scalar intersect(@{ $match_list{$k} }, @{ $match_list{$l} });
                if ( $arrSize > 0 ) {
                        @{ $match_list{$k} } = unique( @{ $match_list{$k} }, @{ $match_list{$l} });
                        undef @{ $match_list{$l} };
                }
        }
}

print "ClusterID\tNum_Of_Species\tB73\tSorghum\tRice\tBrachypodium\tMillet\tBanana\tPalm\tArabidopsis\tAmborella\n";

my $final_count = 1;
for (my $count = 1; $count <= $i; $count++) {
	if ( @{ $match_list{$count} } ) {
		#print "$count\t";
		my $B73_gene_exist = 0;
		#my $W22_gene_exist = 0;
		my $Rice_gene_exist = 0;
		my $Brachy_gene_exist = 0;
		my $Sorghum_gene_exist = 0;
		my $Millet_gene_exist = 0;
		my $Banana_gene_exist = 0;
		my $Palm_gene_exist = 0;
		my $Arabidopsis_gene_exist = 0;
		my $Amborella_gene_exist = 0;
		### find num of species in the ortho group
		foreach my $key ( sort @{ $match_list{$count} } ) {
			#if ($key =~ /Zm/) {
			#	$W22_gene_exist = 1;
			#}
			if ($key =~ /LOC/) {
				$Rice_gene_exist = 1;
			}
			if ($key =~ /Brad/) {
				$Brachy_gene_exist = 1;
			}
			if ($key =~ /GRMZM|AC|AY|EF|AF/) {
				$B73_gene_exist = 1;
			}
			if ($key =~ /Sobic/) {
				$Sorghum_gene_exist = 1;
			}
			if ($key =~ /Si/) {
				$Millet_gene_exist = 1;
			}  		
			if ($key =~ /GSMUA/) {
				$Banana_gene_exist = 1;
			}
			if ($key =~ /p5/) {
				$Palm_gene_exist = 1;
			}
			if ($key =~ /^AT/) {
				$Arabidopsis_gene_exist = 1;
			}
			if ($key =~ /AmTr/) {
				$Amborella_gene_exist = 1;
			}
		}
		my $num_of_species = $B73_gene_exist + $Rice_gene_exist + $Brachy_gene_exist + $Sorghum_gene_exist + $Millet_gene_exist + $Banana_gene_exist + $Palm_gene_exist + $Arabidopsis_gene_exist + $Amborella_gene_exist;
		### print out the conserved AS events in each species, if not exist NA
		### at least two species conserved
		if ($num_of_species > 1) {
			my $final_AS_event_ID = $AS_event_type."_".$final_count."_".$clusterID{$count};
			print "$final_AS_event_ID\t$num_of_species\t";
			### print out B73
			if ($B73_gene_exist == 1) {
				foreach my $event ( sort @{ $match_list{$count} } ) {
					if ($event =~ /GRMZM|AC|AY|EF|AF/) {
						print "$event;"
					}
				}
				print "\t";
			}
			if ($B73_gene_exist == 0) {
				print "NA\t";
			}
			### print out W22
			#if ($W22_gene_exist == 1) {
                	#	foreach my $event ( sort @{ $match_list{$count} } ) {
                    	#		if ($event =~ /Zm/) {
                        #			print "$event;"
                    	#		}
                	#	}
                	#	print "\t";
            		#}
            		#if ($W22_gene_exist == 0) {
                	#	print "NA\t";
           		#}
			### print out Sorghum
			if ($Sorghum_gene_exist == 1) {
                		foreach my $event ( sort @{ $match_list{$count} } ) {
                    			if ($event =~ /Sobic/) {
                        			print "$event;"
                    			}
                		}
                		print "\t";
            		}
            		if ($Sorghum_gene_exist == 0) {
                		print "NA\t";
            		}
			### print out Rice
			if ($Rice_gene_exist == 1) {
                		foreach my $event ( sort @{ $match_list{$count} } ) {
                    			if ($event =~ /LOC/) {
                        			print "$event;"
                    			}
                		}
                		print "\t";
            		}
            		if ($Rice_gene_exist == 0) {
                		print "NA\t";
            		}
			### print out Brachypodium
			if ($Brachy_gene_exist == 1) {
                		foreach my $event ( sort @{ $match_list{$count} } ) {
                    			if ($event =~ /Brad/) {
                        			print "$event;"
                    			}
                		}
                		print "\t";
            		}
            		if ($Brachy_gene_exist == 0) {
                		print "NA\t";
           	 	}
			### print out Millet
			if ($Millet_gene_exist == 1) {
                		foreach my $event ( sort @{ $match_list{$count} } ) {
                    			if ($event =~ /Si/) {
                        			print "$event;"
                    			}
                		}
                		print "\t";
            		}
        		if ($Millet_gene_exist == 0) {
                		print "NA\t";
            		}
			### print out Banana
			if ($Banana_gene_exist == 1) {
                		foreach my $event ( sort @{ $match_list{$count} } ) {
                    			if ($event =~ /GSMUA/) {
                 		       		print "$event;"
                    			}
                		}
                		print "\t";
            		}
        		if ($Banana_gene_exist == 0) {
                		print "NA\t";
            		}
			### print out Palm
			if ($Palm_gene_exist == 1) {
                		foreach my $event ( sort @{ $match_list{$count} } ) {
                    			if ($event =~ /p5/) {
                        			print "$event;"
                    			}
                		}
                		print "\t";
            		}	
            		if ($Palm_gene_exist == 0) {
               	 		print "NA\t";
            		}
           		 ### print out Arabidopsis
            		if ($Arabidopsis_gene_exist == 1) {
                		foreach my $event ( sort @{ $match_list{$count} } ) {
                    			if ($event =~ /^AT/) {
                        			print "$event;"
                    			}
                		}
                		print "\t";
        		}
            		if ($Arabidopsis_gene_exist == 0) {
                		print "NA\t";
            		}
			### print out Amborella
			if ($Amborella_gene_exist == 1) {
                		foreach my $event ( sort @{ $match_list{$count} } ) {
                    			if ($event =~ /AmTr/) {
                        			print "$event;"
                    			}
                		}
                		print "\n";
            		}
            		if ($Amborella_gene_exist == 0) {
                		print "NA\n";
            		}
			$final_count++;			
		}
	}
}

### wublast output tabular format Wublast2.0
#1 	qid 	query sequence identifier
#2 	sid 	subject (database) sequence identifier
#3 	E 	the expectation or E-value
#4 	N 	the number of scores considered jointly in computing E
#5 	Sprime 	the normalized alignment score, expressed in units of bits
#6 	S 	the raw alignment score
#7 	alignlen 	the overall length of the alignment including any gaps
#8 	nident 	the number of identical letter pairs
#9 	npos 	the number of letter pairs contributing a positive score
#10 	nmism 	the number of mismatched letter pairs
#11 	pcident 	percent identity over the alignment length (as a fraction of alignlen)
#12 	pcpos 	percent positive letter pairs over the alignment length (as a fraction of alignlen)
#13 	qgaps 	number of gaps in the query sequence
#14 	qgaplen 	total length of all gaps in the query sequence
#15 	sgaps 	number of gaps in the subject sequence
#16 	sgaplen 	total length of all gaps in the subject sequence
#17 	qframe 	the reading frame in the query sequence (+0 for protein sequences in BLASTP and TBLASTN searches)
#18 	qstart 	the starting coordinate of the alignment in the query sequence
#19 	qend 	the ending coordinate of the alignment in the query sequence
#20 	sframe 	the reading frame in the subject sequence (+0 for protein sequences in BLASTP and BLASTX searches)
#21 	sstart 	the starting coordinate of the alignment in the subject sequence
#22 	send 	the ending coordinate of the alignment in the subject sequence

### header letters of the species
#	LOC	Rice
#	Brad	Brachypodium
#	GRMZM|AC|AY|EF|AF	B73
#	Zm	W22
#	Sobic	Sorghum
#	Si	Millet
#	GSMUA	Banana
#	p5	Palm
#	AT	Arabidopsis
#	AmTr	Amborella

### Prefix for each species in the OrthoFinder
#	B73             B73
#	Rice            OSA
#	Sorghum         SB
#	Millet          SIT
#	Brachypodium    BD
#	Palm            EG
#	Banana          MUSA
#	Arabidopsis     ATH
#	Amborella       AMBO

