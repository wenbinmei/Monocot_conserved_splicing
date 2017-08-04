#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#####################################
# Written by
# May 2015 Wenbin Mei
# Brad's Lab
#####################################

### This scripts take a gtf file generate from cufflinks and filter based on the FPKM and isoform fraction ratio

my ($gtf, $FPKM_value, $frac_value, $help);

GetOptions("input|i=s" => \$gtf, "frac_value|v=s" => \$frac_value, "FPKM_value|f=s" => \$FPKM_value, "help|h" => \$help);

if ($help) {
        print <<EOF;
Usae: perl filter_isoform_fraction_FPKM.pl -i [gtf file] -v [FPKM_value] -f [fraction_value]
                Options
                --input|i: gtf file generate from cufflinks
                --frac_value: fraction value (default 0.03)
                --FPKM_value: FPKM value (default 1)
                --help: print help information
EOF
exit;
}

open (GTF, $gtf) or die "Can not open input file: $gtf\n";
while (<GTF>) {
    chomp;
    my @temp  = split(/\s+/, $_);
    if ($temp[2] eq "transcript") {
       my $Isoform_fraction = $temp[15];
       my $Isoform_FPKM     = $temp[13]; 
          $Isoform_fraction =~ s/["|";]//g;
          $Isoform_FPKM     =~ s/["|";]//g; 
          if ( ($Isoform_fraction >= $frac_value) && ($Isoform_FPKM >= $FPKM_value) ) {
               print "$_\n";
          }
    }
    if ($temp[2] eq "exon") {
       my $fraction_exon = $temp[17];
       my $FPKM_exon     = $temp[15];
          $fraction_exon =~ s/["|";]//g;
          $FPKM_exon     =~ s/["|";]//g;
          if ( ($fraction_exon >= $frac_value) && ($FPKM_exon >= $FPKM_value) ) {
               print "$_\n";
          }
    }
}
close GTF;



          

