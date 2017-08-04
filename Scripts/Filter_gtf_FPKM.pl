#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#####################################
# Written by
# August 2012 Wenbin Mei
# Brad's Lab
#####################################

### This scripts take a gtf file generate from cufflinks and also a FPKM value, and filter out all the transcript with FPKM < value
### also keep the transcripts whether they are full reads support or not

my ($gtf, $fpkm, $support, $help);
GetOptions("Input|i=s" => \$gtf, "FPKM|f=s" => \$fpkm, "Full_Reads_Support|r=s" => \$support, "help|h" => \$help);

if ($help) {
        print <<EOF;
Usae: perl Filter_gtf_FPKM.pl -i [gtf file] -f [value] -r [YES or NO]
                Options
                --Input|i: gtf file generate from cufflinks
                --FPKM|f: FPKM value
                --Full_Reads_Support|r: whether there is RNA-Seq reads suppor the full transcripts (YES or NO)
                --help: print help information
EOF
exit;
}

my %keep_transcript = ( );

open (GTF, $gtf) or die "Can not open input file: $gtf\n";
while (<GTF>) {
    chomp;
    my @temp  = split(/\s+/, $_);
    if ($temp[2] eq "transcript") {

       my $FPKM_transcript = $temp[13];
       $FPKM_transcript =~ s/["|";]//g;

       my $transcript_id = $temp[11];
       $transcript_id =~ s/["|";]//g;

       ### do not require full reads support, only exam the FPKM
       if ($support eq "NO") { 
          if ($FPKM_transcript >= $fpkm) {
               print "$_\n";
               $keep_transcript{$transcript_id}++;
          }
       }

       ### require full reads support
       if ($support eq "YES") {
          my $full_reads_support = $temp[23];
          $full_reads_support =~ s/["|";]//g;
          if ( ($full_reads_support eq "yes") && ($FPKM_transcript >= $fpkm) ) {
                   print "$_\n";
               $keep_transcript{$transcript_id}++;
          }
       }
    }          

    if ($temp[2] eq "exon") {
        my $transcript_id = $temp[11];
        $transcript_id =~ s/["|";]//g;
          
          if ( defined $keep_transcript{$transcript_id} ) {
               print "$_\n";
          }
    }
}
close GTF;



          

