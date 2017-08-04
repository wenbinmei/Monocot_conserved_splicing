#!/usr/bin/perl -w
use strict;

### Jan 26, 2014  Wenbin Mei
### filter out the PASA gtf file single transcript loci

my $usage = "$0 <PASA gtf file>";

my ($File) = @ARGV;
my %link;
#chr1    PASA    transcript      28      3756    .       +       .       gene_id "S1"; transcript_id "asmbl_1";
#chr1    PASA    exon    28      104     .       +       .       gene_id "S1"; transcript_id "asmbl_1";
#chr1    PASA    exon    200     1233    .       +       .       gene_id "S1"; transcript_id "asmbl_1";
#chr1    PASA    exon    1621    3218    .       +       .       gene_id "S1"; transcript_id "asmbl_1";
#chr1    PASA    exon    3301    3756    .       +       .       gene_id "S1"; transcript_id "asmbl_1";

open (FILE, $File) or die "Can not open input file: $File\n";
while (<FILE>) {
        chomp;
        next if /^\s*$/;
        my @line = split (/\s+/, $_);
        if ($line[2] eq "transcript") {
                 my $transID = $line[11];
                 if ($transID =~ /"(\S+)"/) {
                        $transID = $1;
                 }
                 my $geneID = $line[9];
                 if ($geneID =~ /"(\S+)"/) {
                        $geneID = $1;
                 }
                 $link{$geneID}++;
        }
}
close FILE;

open (FILTER, $File) or die "Can not open input file: $File\n";
while (<FILTER>) {
        chomp;
        next if /^\s*$/;
        my @line = split (/\s+/, $_);
        my $geneID = $line[9];
        if ($geneID =~ /"(\S+)"/) {
                $geneID = $1;
        }
        if (  (defined $link{$geneID}) && ($link{$geneID} >= 2)  ) {
                print "$_\n";
        }
}
close FILTER;

