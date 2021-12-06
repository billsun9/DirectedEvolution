#!/usr/bin/perl

use strict;
use warnings;

# evaluate search results based on the number of rounds used for finding the top X-th samples
#
# usage:
# ./this_program Umetsu_ArgPep_combo.py.log.1228.txt Umetsu_ArgPep.val 

# rnd[$i] rounds were used for finding the top $top[$i] candidates 
my @top = (0, 2, 4, 9, 19); 
my @rnd;

my @res;
my @val;

open(IN, $ARGV[0]) or die "cannot read $ARGV[0]\n";
while (my $line=<IN>) {
    chomp $line;
    my @cells = split(/\,/, $line);
    push(@res, $cells[0]);
#    print $cells[0], "\n";
}
close(IN);

open(IN, $ARGV[1]) or die "cannot read $ARGV[1]\n";
while (my $line=<IN>) {
    chomp $line;
    my @cells = split(/\,/, $line);
    push(@val, $cells[0]);
#    print $cells[0], "\n";
}
close(IN);

@val = sort {$b <=> $a} @val;

for (my $i=0; $i<@top; $i++) {
    my $nf = 1;
    for (my $j=0; $j<@res; $j++) {
	if ($res[$j] >= $val[$top[$i]]) {
	    $rnd[$i] = $j;
	    $nf = 0;
	    last
	}
    }
    $nf==1 and $rnd[$i]=@res;
}

print join("\t", @rnd), "\n";

