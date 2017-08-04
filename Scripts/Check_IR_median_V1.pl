#!/usr/bin/perl -w 
use strict;
use warnings;
use Getopt::Long;

### 05/2015 Wenbin Mei
###
### perl Check_IR_median.pl -d depth from the alignment -f file include the location of the intron retention -j splice junction file
### This script calculate three different things:  1) percentage of intron has more than two reads covered; 
###                                                2) median intron coverage;
###                                                3) Intron retention ratio (IRR) ---- median_intron_coverage/(median_intron_coverage + splice junction coverage).
###
### 08/05/2015 Wenbin Mei
### since I made the new scripts Make_PASA_OUTPUT_V2.pl, and the coodinate is exon based and if on the negative strand the coordinate is from large coodinate to small coodinate,  
### so I modify couple place in this scripts and also use the regular splice junction file (start,end and depth)
### The chr name in these three files need to match format !!!!!!
########################################################

my ($base_depth, $coverage, $file, $junction_file, $help);
GetOptions("depth|d=s" => \$base_depth, "file|f=s" => \$file, "junction|j=s" => \$junction_file, "help|h" => \$help);

if ($help) {
        print <<EOF;
Usae: perl Check_RetainedIntron_Coverage.pl -d [depth from alignment] -f [file include the location of intron retention events] -j [splice junction coverage file]
                Options
                --file|f:          file include the location of intron retention event 
                --depth|d:         base pair depth file generated from samtools depth
                --junction|j:      splice junction coverage file
                --help:            print help information
EOF
exit;
}

### depth file example
#chr10   2419    1
#chr10   2420    1
#chr10   2421    1

### intron retention events example
### chr	AS type		strand	clusterID	retained_intron_isoform	start	end	spliced_intron_isoform	start	end

#Chr0	retained_intron	+	82	asmbl_100	NA	NA	asmbl_101	4076804	4076992
#Chr0	retained_intron	-	87	asmbl_112	NA	NA	asmbl_110	1592925	1592847

### check the retained intron isoform only

### splice junction file example (splice junction file is intron coordinate based)
#Chr0    832707  833265  154
#Chr0    833332  834093  12

my %link;
my %IR;
my %junction;

#print "load the depth file!!!\n";
### load the depth file
open (Evidence_File, $base_depth) or die "Can not open input file: $base_depth\n";
while (<Evidence_File>) {
      chomp;
      my @line = split (/\s+/, $_);
      my $Scaffoldname = $line[0];
      my $location = $line[1];
      my $depth = $line[2];
      if ($depth >= 2) {    ### make sure each base has at least two reads support
         $link{$Scaffoldname}{$location}=$depth;
      }
}
close Evidence_File;
#print "finish load the depth file\n";


### load the splice junction file
open (Junction_File, $junction_file) or die "Can not open input file: $junction_file\n";
while (<Junction_File>) {
      chomp;
      next if ($_ =~ /IntronStart/);
      my @line = split (/\s+/, $_);
      my $Scaffoldname      = $line[0];
      my $intron_start      = $line[1];
      my $intron_end        = $line[2];
      my $junction_coverage = $line[3];
      $junction{$Scaffoldname}{$intron_start}{$intron_end} = $junction_coverage;
}  
close Junction_File;

open (FILE, $file) or die "Can not open input file: $file\n";
while (<FILE>) {
      chomp;
      if ($_ =~/AS_type/) {
           print "$_\tPercentage_Coverage\tMedian_Coverage\tIR_Ratio\n";
           next;
      }
      $coverage = 0;
      my @base_coverage = ();
      my $average_coverage = 0;
      my @line = split (/\s+/, $_);      
      my $chr   = $line[0];
      #if ($chr  =~ /scaffold/) { 
      #      $chr   = "evm_27.TU."."$chr";  ### this is only for amborella genome
      #}
      my ($start, $end);
      
      if ($line[1] =~ /retained_intron/) {  ### only check the retained intron here
         
         if ($line[2] =~ /\+/) { ### on the positive strand and convert the location from exon to the intron, make the start to end is from smaller to big which is match the junction file
              $start = $line[8] + 1;
              $end   = $line[9] - 1;
         }
        
         if ($line[2] =~ /\-/) {
              $start = $line[9] + 1;
              $end   = $line[8] - 1;
         }

         for (my $i = $start; $i <= $end; $i++) {
              if ( defined $link{$chr}{$i} ) {
                     $coverage++;                 ##### add up all the base together
                     push (@base_coverage, $link{$chr}{$i});  #### add all the base coverage to a array then can be used to calculat the median
              }
              elsif ( not defined $link{$chr}{$i} ) {
                     push (@base_coverage, '0');
              }
         }

         my $percentage_coverage = $coverage/($end - $start + 1);       
         
         my @vals = sort(@base_coverage);
         #print "SORTED: @vals\n";
         my $sum;
         my $med;
         my $IRR;
         my $junction_coverage = 0;

         ### junction file is intron coordinate based
         if (defined $junction{$chr}{$start}{$end}) {
            $junction_coverage = $junction{$chr}{$start}{$end}; 
         }
 
         ### no splice junction associate with this intron retention event, probably not real
         if (not defined $junction{$chr}{$start}{$end}) {
            next;
         }

         if ( (@vals % 2 == 0) && ($junction_coverage > 0) ){ ### make sure junction coverage > 0
            $sum = $vals[(@vals/2)-1] + $vals[(@vals/2)];
            $med = $sum/2;
            $IRR = $med/($med + $junction_coverage);
            print "$_\t$percentage_coverage\t$med\t$IRR\n";
         }

         if ( (@vals % 2 != 0) && ($junction_coverage > 0) ){ ### make sure junction coverage > 0 
            $med = $vals[@vals/2];
            $IRR = $med/($med + $junction_coverage); 
            print "$_\t$percentage_coverage\t$med\t$IRR\n";
         }  
         
      }
         
}
close FILE;





