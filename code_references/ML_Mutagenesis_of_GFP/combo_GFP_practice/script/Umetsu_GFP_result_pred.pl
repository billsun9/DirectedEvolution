#!/usr/bin/perl

use strict;
use warnings;

# usage:
# ./this_program ST-scale_pred seed npred 

my $RD = "/home/saito/work/160215_aiseq/data/Umetsu_GFP";
my $RR = "/home/saito/work/160215_aiseq/result/Umetsu_GFP/randbase5000";

my $NTrain = 155;

my $FeatName = $ARGV[0]; 
my $SeedFile = $ARGV[1]; 
my $Npred = $ARGV[2];

my $tbl;
my @head;
my %rnk;
my %cnt;

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

    exists($rnk{$name}) and die "multiple entries $name\n";
    $rnk{$name} = 0;
    exists($cnt{$name}) and die "multiple entries $name\n";
    $cnt{$name} = 0;
}
close(IN);

open(IN, "cat $SeedFile | ") or die "cannot read $SeedFile\n";
my @Seed = split("\n", join("", <IN>));
close(IN);

for (my $i=0; $i<@Seed; $i++) {
    my $logfile = "$RR/Umetsu_GFP_${FeatName}/Umetsu_GFP_${FeatName}.log.$Seed[$i].$Npred";

    open(IN, "tail -n10 $logfile | ") or die "cannot read $logfile\n";
    my $tmp = join("", <IN>);
    close(IN);
    $tmp =~ /\[\s*(.*?)\s*\]/s;
    $tmp = $1;
    my @res = split(/\s+/, $tmp);

    for (my $j=0; $j<@res; $j++) {
	my $name = $tbl->[$res[$j]]->[0];
	$rnk{$name} += $j+1;
	$cnt{$name}++;
    }
}

foreach my $name (keys %rnk) {
    $cnt{$name}==0 and next;
    $rnk{$name} /= $cnt{$name};
}

print join(",", @head, "rank", "count"), "\n";
for (my $i=0; $i<@$tbl; $i++) {
    my $name = $tbl->[$i]->[0];
    $cnt{$name}==0 and next;

    print join(",", @{$tbl->[$i]}, $rnk{$name}, $cnt{$name}), "\n";
}
