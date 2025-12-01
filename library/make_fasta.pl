#!/usr/bin/perl
#Eric Morrison
#05172024

use strict;
use warnings;

my $in = $ARGV[0];
open(IN, "$in") || die "Can't open input.\n";

chomp(my @in = <IN>);

if(scalar(@in) == 1){
    $in[0] =~ s/\r|\r\n/\n/g;
    @in = split("\n", $in[0]);
}

my $line_num = 1;

foreach my $line (@in){
    print ">ASV_$line_num\n";
    print $line, "\n";
    $line_num++;
}
