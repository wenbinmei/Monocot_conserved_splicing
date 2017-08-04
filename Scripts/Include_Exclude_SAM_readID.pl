#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

# setup my default

my ($SAM, $help);
my $include_readID = "empty_file";
my $exclude_readID = "empty_file";

GetOptions("input|s=s" => \$SAM, "include|i=s" => \$include_readID, "exclude|e=s" => \$exclude_readID, "help|h" => \$help);

if ($help) {
        print <<EOF;
Usae: perl Include_Exclude_SAM_readID.pl -s [SAM file] -i [readID want to keep] -e [readID want to exclude]
                Options
                --input|s: SAM file
                --include|i: readID want to keep in the alignment
                --exclude|e: readID want to exclude in the alignment
                --help: print help information
EOF
exit;
}

my %keepID_link;
my %removeID_link;

if ($include_readID ne "empty_file") {
   open (Include_FILE, $include_readID) or die "Can not open input file: $include_readID\n";
   while (<Include_FILE>) {
         chomp;
         $keepID_link{$_}++;
   }
   close Include_FILE;
}

if ($exclude_readID ne "empty_file") {
   open (Exclude_FILE, $exclude_readID) or die "Can not open input file: $exclude_readID\n";
   while (<Exclude_FILE>) {
         chomp;
         $removeID_link{$_}++;
   }
   close Exclude_FILE;
}



open (SAM, $SAM) or die "Can not open input file: $SAM\n";
while (<SAM>) {
    chomp;
    ### print the header of the sam file
    if ($_ =~ /^@/) {
       print "$_\n";
       next;
    }
    ### split the sam file
    my ($read_id, $flag, $rname, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq)  = split(/\s+/, $_);

    ### remove unmapped reads
    if ($flag != 4 && $rname ne "*") {
        # only povide ID list to include
        #
        if ($include_readID ne "empty_file") {        
             if (defined $keepID_link{$read_id}) {
               print "$_\n";
             }
        }
        # only provie ID list to remove
        #
        if ( $exclude_readID ne "empty_file" && not defined $removeID_link{$read_id} ) {
         #   if (not defined $removeID_link{$read_id}) {
                print "$_\n";
         #   }
        }
    }
}
close SAM;


















