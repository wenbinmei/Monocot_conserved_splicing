#!/usr/bin/perl
use Bio::SeqIO;
use Getopt::Std;
use strict;
use warnings;

### This scripts filter out the Isoform ID from fasta file and filter the merged GTF file

my %options = ();
my %keep;

#print "Usage: perl filter_GTF_IsoID.pl -f FastaFile -g GTFfile";

getopts("f:g:", \%options);
if(not defined $options{f}){&usage();}
if(not defined $options{g}){&usage();}

my $fasta = $options{f};
my $gtf = $options{g};

my $seqIn  = Bio::SeqIO->new(-file => "$fasta" , '-format' => 'Fasta');

while (my $seqObj = $seqIn->next_seq()){
	my $displayID = $seqObj->display_id;
	#my $newID = $add_name."_".$displayID;
        $keep{$displayID}++
}

#chr1    Cufflinks       exon    1       104     .       +       .       gene_id "XLOC_000001"; transcript_id "TCONS_00000003"; exon_number "1"; oId "CUFF.74.3"; tss_id "TSS1";
#chr1    Cufflinks       exon    200     313     .       +       .       gene_id "XLOC_000001"; transcript_id "TCONS_00000003"; exon_number "2"; oId "CUFF.74.3"; tss_id "TSS1";
#chr1    Cufflinks       exon    422     604     .       +       .       gene_id "XLOC_000001"; transcript_id "TCONS_00000003"; exon_number "3"; oId "CUFF.74.3"; tss_id "TSS1";
#chr1    Cufflinks       exon    1016    1233    .       +       .       gene_id "XLOC_000001"; transcript_id "TCONS_00000003"; exon_number "4"; oId "CUFF.74.3"; tss_id "TSS1";
#chr1    Cufflinks       exon    1621    2082    .       +       .       gene_id "XLOC_000001"; transcript_id "TCONS_00000003"; exon_number "5"; oId "CUFF.74.3"; tss_id "TSS1";
#chr1    Cufflinks       exon    2171    3218    .       +       .       gene_id "XLOC_000001"; transcript_id "TCONS_00000003"; exon_number "6"; oId "CUFF.74.3"; tss_id "TSS1";
#chr1    Cufflinks       exon    3301    3826    .       +       .       gene_id "XLOC_000001"; transcript_id "TCONS_00000003"; exon_number "7"; oId "CUFF.74.3

open (GTF, $gtf) or die "Can not open input file: $gtf\n";
while (<GTF>) {
    chomp;
    my @temp  = split(/\s+/, $_);
    my $transcript_id = $temp[11];
    $transcript_id =~ s/["|";]//g;
    if (defined $keep{$transcript_id}) {
             print "$_\n";
    }
}

