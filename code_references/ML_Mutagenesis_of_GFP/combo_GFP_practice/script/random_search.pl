#!/usr/bin/perl

use strict;
use warnings;

use List::Util;

# perform random search as a baseline for bayesian optimization
# for each round, randomly choose one candidate, and report the best candidate until then 
#
# usage:
# ./this_program Goodman_Science_2013_feat.csv seed 2>&1 | tee Goodman_Science_2013_random_search.pl.log.seed

our $del = ",";

my @val;

open(IN, $ARGV[0]) or die "cannot read $ARGV[0]\n";
<IN>;
while (my $line=<IN>) {
    chomp $line;
    my @cells = split(/$del/, $line);
    push(@val, $cells[$#cells]);
}
close(IN);

srand($ARGV[1] * $ARGV[1] * 1985); 
#srand($ARGV[1] * $ARGV[1]);
my @val_shf = List::Util::shuffle(@val);

my $max = -1000000;
for (my $i=0; $i<@val_shf; $i++) {
    if ($val_shf[$i] > $max) {
	$max = $val_shf[$i];
    }
    print $max, "\n";
}
