#!/usr/bin/perl

use strict;
use warnings;

open(F1, "$ARGV[0]") || die "can't open sam input file: $!";
open(FOUT1,"> $ARGV[1]") || die "cant open output file: $!";

my $chrom=$ARGV[2];
my $left=$ARGV[3];
my $right=$ARGV[4];
my $width=$ARGV[5];
my $dist=$ARGV[6];
my ($first, $second);

while (defined($first = <F1>) && defined($second = <F1>)) {
    my @bed1_line =split (/\t/,$first);
    my @bed2_line = split(/\t/,$second);
    
    my $bed1_bait = $bed1_line[2] eq $chrom && (($bed1_line[3] >= ($left-$width) && $bed1_line[3] <= $left) || ($bed1_line[3] >= $right && $bed1_line[3] <= ($right+$width)));
    my $bed2_bait = $bed2_line[2] eq $chrom && (($bed2_line[3] >= ($left-$width) && $bed2_line[3] <= $left) || ($bed2_line[3] >= $right && $bed2_line[3] <= ($right+$width)));
    
    # trans-interactions only
    my $trans1 = ($bed1_bait && $bed2_line[2] ne $chrom);
    my $trans2 = ( $bed2_bait && $bed1_line[2] ne $chrom);
    
    # cis-interactions
    my $cis1 = ($bed1_bait && $bed2_line[2] eq $chrom && abs($bed1_line[3]-$bed2_line[3])>$dist);
    my $cis2 = ($bed2_bait && $bed1_line[2] eq $chrom && abs($bed1_line[3]-$bed2_line[3])>$dist);
    
    if(($bed1_line[0] eq $bed2_line[0]) && ($bed1_line[4]>0 && $bed2_line[4]>0) && ($trans1 == 1 || $trans2 ==1 || $cis1 == 1 || $cis2 == 1)) {
        print FOUT1 "$bed1_line[2]\t", $bed1_line[3]-10, "\t", $bed1_line[3]+10, "\n";
        print FOUT1 "$bed2_line[2]\t", $bed2_line[3]-10, "\t", $bed2_line[3]+10, "\n";
    }

}

close(FOUT1);
close(F1);

