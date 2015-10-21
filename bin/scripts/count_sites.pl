#!/usr/bin/perl
use strict;
use warnings;
#open(F1,'<',$ARGV[0]) || die "cannot open the input file: $!";

#print "Hello,world\n";

my $size_trans = $ARGV[1];
my $size_cis = $ARGV[2];
my $window_cis = $ARGV[3];

my $bait_chr = $ARGV[4];
my $bait_sta = $ARGV[5];
my $bait_end = $ARGV[6];

my @cis_coord1;
my @cis_coord2;
my $cis_enzyme_num = 0;
my @cis_cut;

my %trans_enzyme_sites;
my %trans_cut_sites;

my %trans_coord1;
my %trans_coord2;
my $trans_enzyme_num = 0;
my %trans_cut;

my $old_chr = "";
my @chroms;

open(F1, '<',$ARGV[0]) || die "cannot open the input file: $!";
open(OUTPUT, '>', $ARGV[7]) || die "cannot open the output file: $!";
open(LOG, '>', 'log.txt') or die "cannot open the log file: $!";
	
while (<F1>) {
	my @file = split;
	my $chr = $file[0];

	## read cis-interactions into array
	if ($chr eq $bait_chr) {
		$cis_coord1[$cis_enzyme_num] = $file[1];
		$cis_coord2[$cis_enzyme_num] = $file[2];
		$cis_cut[$cis_enzyme_num] = $file[3];
		$cis_enzyme_num ++;
	}	

	## read trans-interactions into array
	else {
		
		if ($file[0] ne $old_chr) { # count a new chromosome, reset $trans_enzyme_num = 0
			$old_chr = $file[0];
			$trans_enzyme_num = 0;
			$trans_cut_sites{$file[0]} = 0;
			push @chroms, $old_chr;
		}
				
		$trans_enzyme_sites{$old_chr} ++;
		$trans_cut_sites{$old_chr} = $trans_cut_sites{$old_chr} + $file[3];

		$trans_coord1{$old_chr}[$trans_enzyme_num] = $file[1];
		$trans_coord2{$old_chr}[$trans_enzyme_num] = $file[2];
		$trans_cut{$old_chr}[$trans_enzyme_num] = $file[3];

		$trans_enzyme_num ++;		
	}
}
close(F1);

print "DEBUG:\t@chroms\n";

# calculate z-score for each site in cis
for (my $i=$size_cis/2; $i<$window_cis/2; $i++) {
	my $back_count = 0;
	my $front_count = 0;
	for (my $j=0; $j<=$window_cis; $j++) {
		$back_count = $back_count + $cis_cut[$j];
	}
	for (my $k=$i-$size_cis/2; $k<=$i+$size_cis/2; $k++) {
		$front_count = $front_count + $cis_cut[$k];
	}
	my $p = $back_count/$window_cis;

	if ($p==0) {
		print LOG "The window at $bait_chr\t$cis_coord1[$i-$size_cis/2]\t$cis_coord2[$i+$size_cis/2] has no sites cut. DEBUG: front_count: $front_count\tback_count: $back_count\n";
		next;
	}
	
	my $z_score = ($front_count-$p*$size_cis)/sqrt($p*$size_cis*(1-$p));
	print OUTPUT "$bait_chr\t$cis_coord1[$i-$size_cis/2]\t$cis_coord2[$i+$size_cis/2]\t$z_score\t$i\t$cis_coord1[$i]\t$cis_coord2[$i]\t$front_count\t$back_count\n";
}

for (my $i=$window_cis/2; $i<$cis_enzyme_num-$window_cis/2; $i++) {
	my $back_count = 0;
	my $front_count = 0;
	for (my $j=$i-$window_cis/2; $j<=$i+$window_cis/2; $j++) {
		$back_count = $back_count + $cis_cut[$j];
	}
	for (my $k=$i-$size_cis/2; $k<=$i+$size_cis/2; $k++) {
		$front_count = $front_count + $cis_cut[$k];
	}
	my $p = $back_count/$window_cis;
	
	if ($p==0) {
		print LOG "The window at $bait_chr\t$cis_coord1[$i-$size_cis/2]\t$cis_coord2[$i+$size_cis/2] has no sites cut. DEBUG: front_count: $front_count\tback_count: $back_count\n";
		next;
	}
	
	my $z_score = ($front_count-$p*$size_cis)/sqrt($p*$size_cis*(1-$p));
	print OUTPUT "$bait_chr\t$cis_coord1[$i-$size_cis/2]\t$cis_coord2[$i+$size_cis/2]\t$z_score\t$i\t$cis_coord1[$i]\t$cis_coord2[$i]\t$front_count\t$back_count\n";
}


for (my $i=$cis_enzyme_num-$window_cis/2; $i<$cis_enzyme_num-$size_cis/2; $i++) {
	my $back_count = 0;
	my $front_count = 0;
	for (my $j=$cis_enzyme_num-$window_cis; $j<$cis_enzyme_num; $j++) {
		$back_count = $back_count + $cis_cut[$j];
	}
	for (my $k=$i-$size_cis/2; $k<=$i+$size_cis/2; $k++) {
		$front_count = $front_count + $cis_cut[$k];
	}
	my $p = $back_count/$window_cis;
	
	if ($p==0) {
		print LOG "The window at $bait_chr\t$cis_coord1[$i-$size_cis/2]\t$cis_coord2[$i+$size_cis/2] has no sites cut. DEBUG: front_count: $front_count\tback_count: $back_count\n";
		next;
	}
	
	my $z_score = ($front_count-$p*$size_cis)/sqrt($p*$size_cis*(1-$p));
	print OUTPUT "$bait_chr\t$cis_coord1[$i-$size_cis/2]\t$cis_coord2[$i+$size_cis/2]\t$z_score\t$i\t$cis_coord1[$i]\t$cis_coord2[$i]\t$front_count\t$back_count\n";
}


# calculate z-score for each trans-chromosome
foreach my $chr (@chroms) {
	my $window_trans = $trans_enzyme_sites{$chr};
	my $back_count = $trans_cut_sites{$chr};
	
	my $p = $back_count/$window_trans;
		
	if ($p==0) {
		print LOG "The windows at", $chr, "have no sites cut. DEBUG: back_count: $back_count\n";
		next;
	}
	
	for (my $i=$size_trans/2; $i<$trans_enzyme_sites{$chr}-$size_trans/2; $i++) {
		my $front_count = 0;
		for (my $j=$i-$size_trans/2; $j<=$i+$size_trans/2; $j++) {
			$front_count = $front_count + $trans_cut{$chr}[$j];
		}
		
		my $z_score = ($front_count-$p*$size_trans)/sqrt($p*$size_trans*(1-$p));
		print OUTPUT "$chr\t$trans_coord1{$chr}[$i-$size_trans/2]\t$trans_coord2{$chr}[$i+$size_trans/2]\t$z_score\t$i\t$trans_coord1{$chr}[$i]\t$trans_coord2{$chr}[$i]\t$front_count\t$back_count\n";
	}
}

foreach my $chr (@chroms) {
	print "DEBUG:\t$trans_enzyme_sites{$chr}\t$trans_cut_sites{$chr}\n";
}
