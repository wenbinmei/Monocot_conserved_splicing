#!/usr/bin/perl
use Bio::SeqIO;
use Getopt::Std;
use strict;
use warnings;

my %options = ();

#print "Usage: perl pick_length_fasta.pl -f fastaSeqsFile -l length -o output";

getopts("f:l:o:", \%options);
if(not defined $options{f}){&usage();}
if(not defined $options{l}){&usage();}
if(not defined $options{o}){&usage();}

my $seqsFile = $options{f};
my $length   = $options{l};
my $outFile   = $options{o};
my $seqlen;

my $seqIn  = Bio::SeqIO->new(-file => "$seqsFile" , '-format' => 'Fasta');
my $seqOut = Bio::SeqIO->new(-file => ">$outFile" , '-format' => 'Fasta');

while (my $seqObj = $seqIn->next_seq()){
	my $displayID = $seqObj->display_id;
	$seqlen = $seqObj->length;
        if ($seqlen >= $length) {
             $seqOut->write_seq($seqObj);
        }
}

