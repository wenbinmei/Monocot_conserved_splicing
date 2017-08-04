#!/usr/bin/perl -w
use strict;

### this script take an input file GRMZMONYMs_SYNTELOGs_FrJamesSchnable_Sep2015.csv and 
### output the file with 5 species only maize, sorghum, rice, brachypodium, foxtail millet
### select synteny with gene exist in all five species, either maize1, maize2 or both

#maize1,maize2,sorghum,setaria,panicA,panicB,rice,brachy,GEvo Link
#GRMZM2G007948,No Gene,Sobic.001G000100,No Gene,No Gene,No Gene,No Gene,No Gene,http://genomevolution.org/CoGe/GEvo.pl?accn1=Sobic.001G000100;accn2=GRMZM2G007948;accn3=Sevir.9G001200.1.v1;accn4=Oropetium_20150105_23638A;num_seqs=4;autogo=1
#GRMZM2G007675,No Gene,Sobic.001G000200,Si035735m,No Gene,No Gene,No Gene,Bradi1g00450,http://genomevolution.org/CoGe/GEvo.pl?accn1=Sobic.001G000200;accn2=GRMZM2G007675;accn3=Si035735m;accn4=Sevir.9G001300.1.v1;accn5=Oropetium_20150105_23637A;accn6=Bradi1g00450;num_seqs=6;autogo=1

my $usage = "$0 <GRMZMONYMs_SYNTELOGs_FrJamesSchnable_Sep2015.csv>";
die $usage unless $#ARGV == 0;
my ($file) = @ARGV;
my %link;

open (FILE, $file) or die "Can not open input file: $file\n";
while (<FILE>) {
          chomp;
          my $count = 0;
          my ($maize1, $maize2, $sorghum, $setaria, $panicA, $panicB, $rice, $brachy, $GEvo) = split(/\,/, $_);
          if ($maize1 !~ /No/) {
          	$count++;
          }
          if ($maize2 !~ /No/) {
          	$count++;
          }
          if ($sorghum !~ /No/) {
          	$count++;
          }
          if ($setaria !~ /No/) {
          	$count++;
          }
          if ($rice !~ /No/) {
          	$count++;
          }
          if ($brachy !~ /No/) {
          	$count++;
          }
          ### maize1 missing
          if ($maize1 =~ /No/ && $count >= 5) {
          	#print $maize1.",".$maize2.",".$sorghum.",".$setaria.",".$rice.",".$brachy."\n";
          	$link{$maize2}++;
          	$link{$sorghum}++;
          	$link{$setaria}++;
          	$link{$rice}++;
          	$link{$brachy}++;
          }
          ### maize2 missing
          if ($maize2 =~ /No/ && $count >= 5) {
          	#print $maize1.",".$maize2.",".$sorghum.",".$setaria.",".$rice.",".$brachy."\n";
          	$link{$maize1}++;
          	$link{$sorghum}++;
          	$link{$setaria}++;
          	$link{$rice}++;
          	$link{$brachy}++;
          }
          ### both maize1 maize2 exist
          if ($count == 6) {
          	#print $maize1.",".$maize2.",".$sorghum.",".$setaria.",".$rice.",".$brachy."\n";
          	$link{$maize1}++;
          	$link{$maize2}++;
          	$link{$sorghum}++;
          	$link{$setaria}++;
          	$link{$rice}++;
          	$link{$brachy}++;
          }
}
close FILE;

open (OUTPUT, $file) or die "Can not open input file: $file\n";
while (<OUTPUT>) {
          chomp;
          my $count = 0;
          if ($_ =~ /maize/) {
          	print "maize1,maize2,sorghum,setaria,rice,brachy\n";
          	next;
          }
          my ($maize1, $maize2, $sorghum, $setaria, $panicA, $panicB, $rice, $brachy, $GEvo) = split(/\,/, $_);
          if ($maize1 !~ /No/) {
          	$count++;
          }
          if ($maize2 !~ /No/) {
          	$count++;
          }
          if ($sorghum !~ /No/) {
          	$count++;
          }
          if ($setaria !~ /No/) {
          	$count++;
          }
          if ($rice !~ /No/) {
          	$count++;
          }
          if ($brachy !~ /No/) {
          	$count++;
          }
          ### maize1 missing
          if ($maize1 =~ /No/ && $count >= 5) {
          	if ( ($link{$maize2} == 1) && ($link{$sorghum} == 1) && ($link{$setaria} == 1) && ($link{$rice} == 1) && ($link{$brachy} == 1) ) {
				print $maize1.",".$maize2.",".$sorghum.",".$setaria.",".$rice.",".$brachy."\n";
		}
          }
          ### maize2 missing
          if ($maize2 =~ /No/ && $count >= 5) {
          	if ( ($link{$maize1} == 1) && ($link{$sorghum} == 1) && ($link{$setaria} == 1) && ($link{$rice} == 1) && ($link{$brachy} == 1) ) {	
				print $maize1.",".$maize2.",".$sorghum.",".$setaria.",".$rice.",".$brachy."\n";
          	}
          }
          ### both maize1 maize2 exist
          if ($count == 6) {
		if ( ($link{$maize1} == 1) && ($link{$maize2} == 1) && ($link{$sorghum} == 1) && ($link{$setaria} == 1) && ($link{$rice} == 1) && ($link{$brachy} == 1) ) {
          			print $maize1.",".$maize2.",".$sorghum.",".$setaria.",".$rice.",".$brachy."\n";
          	}
          }
}
close OUTPUT;

