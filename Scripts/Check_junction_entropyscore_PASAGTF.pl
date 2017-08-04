#!/usr/bin/perl -w 
use strict;
use warnings;
use Getopt::Long;


### This script identify the list of transcript ID which does not have enough suppor for each junction in the transcript
### perl Check_junctionSupport_CufflinksGTF.pl -g NS19_transcripts.gtf -i NS19_transcripts.intronSize -s juncs.all(produced from spanki) -c 2

my ($GTF, $junction_GTF, $junction_SAM, $coverage, $help);
GetOptions("GTF|g=s" => \$GTF, "Junction_GTF|i=s" => \$junction_GTF, "Junction_SAM|s=s" => \$junction_SAM, "coverage|c=s" => \$coverage, "help|h" => \$help);

if ($help) {
        print <<EOF;
Usae: perl Check_junction_entropyscore_PASAGTF.pl -g [gtf file generated from cufflinks] -i [Junction file generate from GTF] -s [Junction generate from spanki] -c [coverage]
                Options
                --GTF|g: transcript.gtf file generated from PASA
                --Junction_GTF|i: junction file generated from parse_intron.pl based on PASA transcript.gtf
                --Junction_SAM|s: junction file generated from spanki from alignment
                --coverage|c: entropy score for each junction
                --help: print help information
EOF
exit;
}


# junction file from spanki
# juncid  dinucleotide    intron_size     annostatus      gmcode  regcode geneassign      geneassignL     geneassignR     unfilt_cov      cov     normcov offsets entropy hamming3        hamming5        MAXmmes MAXminanc
# chr10:100012751_100012933:+     AA..AG  183     un      1t.0i.1s        se      GRMZM2G154667   GRMZM2G154667   GRMZM2G154667   1       1       0.00110367097931        1       0.000000        8       8       17      47
# chr10:100014945_100018924:-     CA..TT  3980    un      -       -       ambiguous       GRMZM2G455420   none    1       1       0.00110367097931        1       0.000000        8       9       2       16


# Location for each junction extract from transcript.gtf Cufflinks
### Intron_Name     Chr     Start   End     Size
### GRMZM2G093399.GRMZM2G093399_T01.intron1/3       chr1    136635  136718  83

my %link;
my %splice_site_number;
my $transID = " ";
my %Filter;

open (Evidence_File, $junction_SAM) or die "Can not open input file: $junction_SAM\n";

while (<Evidence_File>) {
      chomp;
      next if (/juncid/);
      my @line = split (/\s+/, $_);
      #my @line = split (/\t/, $_);
      my $id = $line[0];
      my @info = split (/\:/, $id);
      my $Scaffoldname = $info[0];
      my ($start, $end) = split (/\_/, $info[1]);
      my $entropy = $line[13];
      #print "$Scaffoldname\t$start\t$end\t$entropy\n";
      if ($entropy >= $coverage) {
             $link{$Scaffoldname}{$start}{$end}=$entropy;
      }
}
close Evidence_File;


open (Predict_File, $junction_GTF) or die "Can not open input file: $junction_GTF\n";
while (<Predict_File>) {
      next if (/^Intron_Name/);
      chomp;
      my @line = split (/\s+/, $_);
      if ($line[0] =~ /(asmbl_\d+).intron/) {
                 $transID = $1;
      }
      if ($line[0] =~ /(align_id:\S+).intron/) {
                 $transID = $1;
      }
  #    if ($line[0] =~ /(GRMZM2G\d+_T\d+).intron/) {
  #              $transID = $1;
  #    }
  #    elsif ($line[0] =~ /(AC\d+\.\d+_FGT\d+).intron/) {
  #              $transID = $1;
  #    }
  #    elsif ($line[0] =~ /(CUFF\.\d+\.\d+).intron/) {
  #              $transID = $1;
         #       print "$transID\n";
  #    }
      my $Scaffoldname = $line[1];
      my $start = $line[2];
      my $end = $line[3];      
      if (not defined $link{$Scaffoldname}{$start}{$end}) {
            print "$transID\t$Scaffoldname\t$start\t$end\n";
            $Filter{$transID}++;
      }
}
close Predict_File;

my $OUTPUT = $GTF."_afterfilter_entropyscore"."$coverage".".gtf";

open (OUTPUT, ">$OUTPUT");
open (Transcript, $GTF) or die "Can not open input file: $GTF\n";
while (<Transcript>) {
      chomp;
      next if /^\s*$/;
      my @line = split (/\s+/, $_);
      my $transID = $line[11];  
      if ($transID =~ /"(\S+)"/) {
            $transID = $1;
      }
      if (defined $Filter{$transID}) {
            next;
      }    
      print OUTPUT "$_\n";
}
close Transcript;











