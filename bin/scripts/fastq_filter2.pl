#!/usr/bin/perl
use List::Util qw(sum);

use strict;
use warnings;


my @string;
my @line;
my $count = 0;
my $mean1 = 0;
my $mean2 = 0;
my @score1;
my @score2;

while (<>) {
    chomp;
    $count++;
    if ($count %4 != 0){$line[$count]=$_;}
    else {
        @string = split/\t/,$_;
        @score1 = split//,$string[0];
        @score2 = split//,$string[1];
        @score1 = map(ord,@score1);
        @score2 = map(ord,@score2);
        $mean1 = sum(@score1)/@score1;
        $mean2 = sum(@score2)/@score2;
        if($mean1 >= 53 && $mean2 >= 53) {
            print $line[1],"\n",$line[2],"\n",$line[3],"\n",$_,"\n";
        }
        $count = 0;
    }

}
