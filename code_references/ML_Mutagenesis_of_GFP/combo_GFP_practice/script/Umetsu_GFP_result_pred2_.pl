#!/usr/bin/perl

use strict;
use warnings;

# usage:
# ./this_program ST-scale_pred seed npred 

my $RD = "/home/saito/work/160215_aiseq/data/Umetsu_GFP";
my $RR = "/home/saito/work/160215_aiseq/result/Umetsu_GFP/randbase5000";

my $NTrain = 154;

my $FeatName = $ARGV[0]; 
my $SeedFile = $ARGV[1]; 
my $Npred = $ARGV[2];

my $tbl;
my @head;
my %ff;

my $FeatFile = "$RD/Umetsu_GFP_${FeatName}.csv";
my $n = 0;
open(IN, $FeatFile) or die "cannot read $FeatFile\n";
while (my $line=<IN>) {
    chomp $line;

    if (@head == 0) {
	@head = split(/\,/, $line);
	next;
    }

    $n++ < $NTrain and next;
#    print STDERR "load $line\n";

    my @cells = split(/\,/, $line);
    my $name = shift(@cells);
    my $score = pop(@cells);
    my @feat = @cells;
    push(@$tbl, [$name, @feat, $score]);

    exists($ff{$name}) and die "multiple entries $name\n";
    $ff{$name} = 0;
}
close(IN);

open(IN, "cat $SeedFile | ") or die "cannot read $SeedFile\n";
my @Seed = split("\n", join("", <IN>));
close(IN);

for (my $i=0; $i<@Seed; $i++) {
    my $logfile = "$RR/Umetsu_GFP_${FeatName}/Umetsu_GFP_${FeatName}.log.$Seed[$i].$Npred";

    open(IN, $logfile) or die "cannot read $logfile\n";
    while (my $line=<IN>) {
	$line =~ /^ff\:$/ and last;
    }
    while (my $line=<IN>) {
	$line =~ /^end$/ and last;
	chomp $line;
	my ($j, $f) = split(/\t/, $line);
	my $name = $tbl->[$j]->[0];
	$ff{$name} += $f;
    }
    close(IN);
}

foreach my $name (keys %ff) {
    $ff{$name} /= @Seed;
}

print join(",", @head, "ff"), "\n";
for (my $i=0; $i<@$tbl; $i++) {
    my $name = $tbl->[$i]->[0];
    
    print join(",", @{$tbl->[$i]}, $ff{$name}), "\n";
}
