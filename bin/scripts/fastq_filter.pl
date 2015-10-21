#!/usr/bin/perl
use List::Util qw(sum);

use strict;
use warnings;


my @string;
my @line;
my $count = 0;
my $mean = 0;
while (<>) {
    chomp;
    $count++;
    if ($count %4 ne 0){
        $line[$count]=$_;
    }
    else {
        @string = split//,$_;
        @string = map(ord,@string);
        $mean = sum(@string)/@string;
        if ($mean >= 53) {
            print $line[1],"\n",$line[2],"\n",$line[3],"\n",$_,"\n";
        }
        $count = 0;
    }

}
