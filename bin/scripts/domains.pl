#!/usr/bin/perl
use strict;
use warnings;
open(F1,'<',$ARGV[0]) || die "cannot open the input file: $!";
open(OUTPUT, '>', $ARGV[1]) or die "cannot open the input file: $!";

my $old_chr = "";
my $old_zScore = 9999999999999999;

my $left_coord = 9999999999999999;
my $right_coord = 0;

while (<F1>) {
	my @file = split;
        
	my $chr = $file[0];
        my $zScore = $file[3];
        
        
        if ($old_chr eq "") {#start a new chromosome, from scratch           
            $old_chr = $file[0];
            $old_zScore = $file[3];
            $left_coord = $file[1];
            $right_coord = $file[2];    
        }
        
        if (($file[0] ne $old_chr) && ($old_chr ne "") ) {#start a new chromosome, before further actions, first print 
            print OUTPUT "$old_chr\t$left_coord\t$right_coord\t$old_zScore\n";
            
            $old_chr = $file[0];
            $old_zScore = $file[3];
            $left_coord = $file[1];
            $right_coord = $file[2];    
        }
        if (($file[0] eq $old_chr) && ($file[3] == $old_zScore)) {#same chromosome, same z_score
            $right_coord = $file[2];
        }
        if (($file[0] eq $old_chr) && ($file[3] != $old_zScore)) {#same chromosome, a new z_score
            print OUTPUT "$old_chr\t$left_coord\t$right_coord\t$old_zScore\n";
            $left_coord = $file[1];
            $right_coord = $file[2];
            $old_zScore = $file[3];

        }
}
print OUTPUT "$old_chr\t$left_coord\t$right_coord\t$old_zScore\n";
