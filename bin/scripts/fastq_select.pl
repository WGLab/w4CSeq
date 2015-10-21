#!/usr/bin/perl
use warnings;
use strict;


if($ARGV[4] eq "yes") {
    open(F1, "$ARGV[0]") || die "can't open input fastq file: $!";
}
else {
    open(F1, "gunzip -c $ARGV[0] |") || die "can't open input fastq file: $!";
}

open(FOUT1, '> ', $ARGV[1]) || die "can't open output file: $!";

my $align=$ARGV[2];
my $enzyme=$ARGV[3];
my $start=length($align)-length($enzyme);
my $line1;
my $line2;
my $line3;
my $line4;
my $c=0;

while(my $line =<F1>) {
    $c=$c+1;

    if ($c == 1) {$line1=$line;}
    if ($c == 2) {$line2=$line;}
    if ($c == 3) {$line3=$line;}
    if ($c == 4) {
        $line4=$line;
        if ($line2 =~ /^$align/) {
            print FOUT1 $line1, substr($line2,$start,length($line2)-$start+1),$line3,substr($line4,$start,length($line2)-$start+1);
        }
	$c=0;
    }
}
close(FOUT1);
close(F1);

