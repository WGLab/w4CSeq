#!/var/www/html/w4cseq/bin/localPerl-5-10.1/bin/perl

use strict;
use warnings;

use lib '/var/www/html/w4cseq/bin/localPerl-5-10.1/lib/perl5/site_perl/5.10.1/x86_64-linux/';
use Math::CDF qw(:all);

open(F1,"$ARGV[0]") || die "can't open input file: $!";
open(FOUT1,"> $ARGV[1]") || die "cant open output file: $!";
open(LOG, '>', 'log.txt') or die "cannot open the log file: $!";


my $build = $ARGV[2];
my $size_trans = $ARGV[3];
my $size_cis = $ARGV[4];
my $window_cis = $ARGV[5];
my $bait_chr = $ARGV[6];

my %length_chr;

if ($build eq "mm10"){
	$length_chr{"chr1"}= 195471971;
	$length_chr{"chr2"}= 182113224;
	$length_chr{"chr3"}= 160039680;
	$length_chr{"chr4"}= 156508116;
	$length_chr{"chr5"}= 151834684;
	$length_chr{"chr6"}= 149736546;
	$length_chr{"chr7"}= 145441459;
	$length_chr{"chr8"}= 129401213;
	$length_chr{"chr9"}= 124595110;
	$length_chr{"chr10"}= 130694993;
	$length_chr{"chr11"}= 122082543;
	$length_chr{"chr12"}= 120129022;
	$length_chr{"chr13"}= 120421639;
	$length_chr{"chr14"}= 124902244;
	$length_chr{"chr15"}= 104043685;
	$length_chr{"chr16"}= 98207768;
	$length_chr{"chr17"}= 94987271;
	$length_chr{"chr18"}= 90702639;
	$length_chr{"chr19"}= 61431566;
	$length_chr{"chrX"}= 171031299;
}
if ($build eq "mm9") {
	$length_chr{"chr1"}= 197195432;
	$length_chr{"chr2"}= 181748087;
	$length_chr{"chr3"}= 159599783;
	$length_chr{"chr4"}= 155630120;
	$length_chr{"chr5"}= 152537259;
	$length_chr{"chr6"}= 149517037;
	$length_chr{"chr7"}= 152524553;
	$length_chr{"chr8"}= 131738871;
	$length_chr{"chr9"}= 124076172;
	$length_chr{"chr10"}= 129993255;
	$length_chr{"chr11"}= 121843856;
	$length_chr{"chr12"}= 121257530;
	$length_chr{"chr13"}= 120284312;
	$length_chr{"chr14"}= 125194864;
	$length_chr{"chr15"}= 103494974;
	$length_chr{"chr16"}= 98319150;
	$length_chr{"chr17"}= 95272651;
	$length_chr{"chr18"}= 90772031;
	$length_chr{"chr19"}= 61342430;
	$length_chr{"chrX"}= 166650296;
}

if ($build eq "hg19"){
	$length_chr{"chr1"}= 249250621;
	$length_chr{"chr2"}= 243199373;
	$length_chr{"chr3"}= 198022430;
	$length_chr{"chr4"}= 191154276;
	$length_chr{"chr5"}= 180915260;
	$length_chr{"chr6"}= 171115067;
	$length_chr{"chr7"}= 159138663;
	$length_chr{"chr8"}= 146364022;
	$length_chr{"chr9"}= 141213431;
	$length_chr{"chr10"}= 135534747;
	$length_chr{"chr11"}= 135006516;
	$length_chr{"chr12"}= 133851895;
	$length_chr{"chr13"}= 115169878;
	$length_chr{"chr14"}= 107349540;
	$length_chr{"chr15"}= 102531392;
	$length_chr{"chr16"}= 90354753;
	$length_chr{"chr17"}= 81195210;
	$length_chr{"chr18"}= 78077248;
	$length_chr{"chr19"}= 59128983;
	$length_chr{"chr20"}= 63025520;
	$length_chr{"chr21"}= 48129895;
	$length_chr{"chr22"}= 51304566;
	$length_chr{"chrX"}= 155270560; 
}

if ($build eq "hg18") {
	$length_chr{"chr1"}= 247249719;
	$length_chr{"chr2"}= 242951149;
	$length_chr{"chr3"}= 199501827;
	$length_chr{"chr4"}= 191273063;
	$length_chr{"chr5"}= 180857866;
	$length_chr{"chr6"}= 170899992;
	$length_chr{"chr7"}= 158821424;
	$length_chr{"chr8"}= 146274826;
	$length_chr{"chr9"}= 140273252;
	$length_chr{"chr10"}= 135374737;
	$length_chr{"chr11"}= 134452384;
	$length_chr{"chr12"}= 132349534;
	$length_chr{"chr13"}= 114142980;
	$length_chr{"chr14"}= 106368585;
	$length_chr{"chr15"}= 100338915;
	$length_chr{"chr16"}= 88827254;
	$length_chr{"chr17"}= 78774742;
	$length_chr{"chr18"}= 76117153;
	$length_chr{"chr19"}= 63811651;
	$length_chr{"chr20"}= 62435964;
	$length_chr{"chr21"}= 46944323;
	$length_chr{"chr22"}= 49691432;
	$length_chr{"chrX"}= 154913754; 
}


my %trans_sites;
my %trans_coord;

my $old_chr = "";
my $cut = 0;

my @cis_coord;
my $cis_num=0;
my $trans_num;
my @chroms;

while (<F1>) {
	my @file = split;
	my $chr = $file[0];
	
	if ($chr eq $bait_chr) {
		$cis_coord[$cis_num] = $file[1];
		$cis_num ++;
	}
	
	else {
		if ($file[0] ne $old_chr) { # count a new chromosome, reset $trans_num = 0
			$old_chr = $file[0];
			$trans_num = 0;
			if ($old_chr ne "chrM") {
				push @chroms, $old_chr;
			}
			
		}
		
		$trans_sites{$old_chr} ++;
		$trans_coord{$old_chr}[$trans_num] = $file[1];
		$trans_num ++;
	}
	
	
}

##DEBUG
foreach my $chr (@chroms) {
	print "DEBUG: ", $chr,"\t", $trans_sites{$chr},"\n";
	print "DEBUG: ", $trans_coord{$chr}[$trans_sites{$chr}-1], "\n";
}	

# calculate p-value for each trans-chromosome;
foreach my $chr (@chroms) {
	my $p = $trans_sites{$chr}/$length_chr{$chr};
		
	if ($p==0) {
		print LOG "The windows at $chr have no sites cut. DEBUG: back_count: $trans_sites{$chr}\n";
		next;
	}
	
	my $mu = $p*$size_trans;
	
	for (my $i=0; $i<$trans_sites{$chr}; $i++) {
		my $count = 0;
		
		#left end
		if ($trans_coord{$chr}[$i] < $size_trans/2) {
			for (my $j=0; $j<$trans_sites{$chr}; $j++) {
				if (($trans_coord{$chr}[$j] <= $size_trans/2) || ($trans_coord{$chr}[$j]-$trans_coord{$chr}[$i] <= $size_trans/2+$size_trans/2-$trans_coord{$chr}[0])) {
					$count ++;
				}
			}
		}
		
		#right end
		elsif ($length_chr{$chr}-$trans_coord{$chr}[$i] < $size_trans/2) {
			for (my $j=0; $j<$trans_sites{$chr}; $j++) {
				if (($length_chr{$chr}-$trans_coord{$chr}[$j] <= $size_trans/2) || ($trans_coord{$chr}[$i]-$trans_coord{$chr}[$j] <= $size_trans/2 + $size_trans/2-($length_chr{$chr}-$trans_coord{$chr}[$trans_sites{$chr}-1]))) {
					$count ++;
				}
				
			}
		}
		
		#middle part
		else {
			for (my $j=0; $j<$trans_sites{$chr}; $j++) {
				if (abs($trans_coord{$chr}[$i]-$trans_coord{$chr}[$j]) <= $size_trans/2) {
					$count ++;
				}
				
			}
		}
		if ($count == 0) {
			print FOUT1 $chr, "\t", $trans_coord{$chr}[$i], "\t", $trans_coord{$chr}[$i]+20, "\t1\n";
		}
		else {
			print FOUT1 $chr, "\t", $trans_coord{$chr}[$i], "\t", $trans_coord{$chr}[$i]+20, "\t", 1-pbinom($count-1, $size_trans, $p), "\n";
		}
	}
	
}


# calculate p-value for cis-chromosome
for (my $i=0; $i< scalar @cis_coord; $i++) {
	my $back_count = 0;
	my $front_count = 0;
	
	# left end
	if ($cis_coord[$i] <= $window_cis/2) {
		for (my $j=0; $j< scalar @cis_coord; $j++) {
			if ($cis_coord[$j] <= $window_cis) {
				$back_count ++;
			}	
		}
	}
	
	
	# right end
	elsif ($length_chr{$bait_chr}-$cis_coord[$i] <= $window_cis/2) {
		for (my $j=0; $j< scalar @cis_coord; $j++) {
			if ($cis_coord[$j] >= $length_chr{$bait_chr}-$window_cis) {
				$back_count ++;
			}	
		}
	}
	
	else {
		for (my $j=0; $j< scalar @cis_coord; $j++) {
			if (abs($cis_coord[$i]-$cis_coord[$j]) <= $window_cis/2) {
				$back_count ++;
			}	
		}
	}
	
	if ($back_count == 0) {
		print LOG "The window at $bait_chr\t$cis_coord[$i] has no sites cut. DEBUG: back_count: $back_count\n";
		next;
	}
	
	
	for (my $k=0; $k< scalar @cis_coord; $k++) {
		if (abs($cis_coord[$i]-$cis_coord[$k]) <= $size_cis/2) {
			$front_count ++;
		}
		
	}
	my $p = $back_count/$window_cis;
	
	if ($front_count == 0) {
		print FOUT1 $bait_chr, "\t", $cis_coord[$i], "\t", $cis_coord[$i]+20, "\t1\n";
	}
	else {
		print FOUT1 $bait_chr, "\t", $cis_coord[$i], "\t", $cis_coord[$i]+20, "\t", 1-pbinom($front_count-1, $size_cis, $p), "\n";
	}
}

close (F1);
close (FOUT1);
close (LOG);

