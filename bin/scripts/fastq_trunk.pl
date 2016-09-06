#!/usr/bin/perl
use strict;
use warnings;

if ($ARGV[2] eq "no") {
    open(F1, "gunzip -c $ARGV[0] |") || die "can't open input fastq file: $!";
}

else {
    open(F1, "$ARGV[0]") || die "can't open input fastq file: $!";
}

open(FOUT1, "> $ARGV[1]") || die "can't open output file: $!";

my $line1;
my $line2;
my $line3;
my $line4;
my $c=0;
while(my $line =<F1>) {
    $c=$c+1;

    if ($c == 1) {
	$line1 = $line;
	if ($line1 =~ /\A@(.*\s)\w.*/) { #to be compatible with SRA file format
            $line1 =~ s/$1//;
        }
    }
    if ($c == 2) {$line2 = $line;}
    if ($c == 3) {$line3 = $line;}
    if ($c == 4) {
        $line4 = $line; 
        print FOUT1 $line1, substr($line2,0,20),"\n", $line3, substr($line4,0,20),"\n";
	$c=0;
    }
}
close(FOUT1);
close(F1);

