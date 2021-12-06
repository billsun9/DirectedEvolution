#!/usr/bin/perl

use strict;
use warnings;

# usage:
# ./this_program ST-scale seed 

my $RD = "/home/saito/work/160215_aiseq/data/Umetsu_GFP";
my $RR = "/home/saito/work/160215_aiseq/result/Umetsu_GFP/randbase0";

my $FeatName = $ARGV[0]; 
my $SeedFile = $ARGV[1]; 

my @score;
my $n = 100;
for (my $i=0; $i<$n; $i++) {
    $score[$i] = 0;
}

open(IN, "cat $SeedFile | ") or die "cannot read $SeedFile\n";
my @Seed = split("\n", join("", <IN>));
close(IN);

for (my $i=0; $i<@Seed; $i++) {
    my $logfile = "$RR/Umetsu_GFP_${FeatName}/Umetsu_GFP_${FeatName}.log.$Seed[$i]";

    my $j = 0;
    open(IN, "cat $logfile | grep current | ") or die "cannot read $logfile\n";
    while (my $line=<IN>) {
	if ($line =~ /\=\s(\S+)/) {
	    $score[$j++] += $1;
	}
    }
    close(IN);
    $j==$n or die "wrong log file $logfile\n";
}

for (my $i=0; $i<$n; $i++) {
    $score[$i] /= @Seed;
    print $score[$i], "\n";
}
