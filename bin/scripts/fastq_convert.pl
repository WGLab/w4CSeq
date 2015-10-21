#!/usr/bin/perl
use List::Util 'max';

use strict;
use warnings;


my @string;
my $count = 0;
my $max = 0;
my $flag = 0;
while (<>) {
    chomp;
    if ($count++ % 4 == 3) {
        @string = split//,$_;
        @string = map(ord,@string);
        $max = max(@string);
        if($max > 74) {$flag=1;}
        if($flag == 1) {tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/;}
     }
    print "$_\n";
}
