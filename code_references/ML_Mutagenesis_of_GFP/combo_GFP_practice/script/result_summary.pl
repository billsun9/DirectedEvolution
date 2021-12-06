#!/usr/bin/perl

use strict;
use warnings;

#our $NR=7;
#our $FeatFile = "/home/saito/work/160215_aiseq/data/Umetsu_ArgPep/Umetsu_ArgPep_1mer_2mer_3mer_pred_r$NR.csv";
#our $Prefix = "/home/saito/work/160215_aiseq/result/Umetsu_ArgPep/randbase0/Umetsu_ArgPep_1mer_2mer_3mer_pred_r$NR/Umetsu_ArgPep_1mer_2mer_3mer_pred_r$NR.log";
our $FeatFile = "/home/saito/work/160215_aiseq/data/Umetsu_ArgPep/Umetsu_ArgPep_1mer_2mer_3mer_pred.csv";
our $Prefix = "/home/saito/work/160215_aiseq/result/Umetsu_ArgPep/randbase0/Umetsu_ArgPep_1mer_2mer_3mer_pred/Umetsu_ArgPep_1mer_2mer_3mer_pred.log";

#our $Suffix = 20;
our $Suffix = 100;
our $NTrain = 90;

our @Seed = (
0,
);

my $tbl;
my @head;
my %cnt;

my $n = 0;
open(IN, $FeatFile) or die "cannot read $FeatFile\n";
while (my $line=<IN>) {
    chomp $line;

    if (@head == 0) {
	@head = split(/\,/, $line);
	next;
    }

    $n++ < $NTrain and next;
    print STDERR "load $line\n";

    my @cells = split(/\,/, $line);
    my $name = shift(@cells);
    my $score = pop(@cells);
    my @feat = @cells;
    push(@$tbl, [$name, @feat, $score]);

    exists($cnt{$name}) and die "multiple entries $name\n";
    $cnt{$name} = 0;
}
close(IN);


for (my $i=0; $i<@Seed; $i++) {
    my $logfile = join(".", ($Prefix, $Seed[$i], $Suffix));
    my $resfile = join(".", ($Prefix, $Seed[$i], $Suffix, "csv"));

    open(IN, "tail -n10 $logfile | ") or die "cannot read $logfile\n";
    my $tmp = join("", <IN>);
    close(IN);
    $tmp =~ /\[\s*(.*?)\s*\]/s;
    $tmp = $1;
    my @res = split(/\s+/, $tmp);

    for (my $j=0; $j<@res; $j++) {
	my $name = $tbl->[$res[$j]]->[0];
	$cnt{$name}++;
    }
}

for (my $i=0; $i<@Seed; $i++) {
    my $logfile = join(".", ($Prefix, $Seed[$i], $Suffix));
    my $resfile = join(".", ($Prefix, $Seed[$i], $Suffix, "txt"));

    open(IN, "tail -n10 $logfile | ") or die "cannot read $logfile\n";
    my $tmp = join("", <IN>);
    close(IN);
    $tmp =~ /\[\s*(.*?)\s*\]/s;
    $tmp = $1;
    my @res = split(/\s+/, $tmp);

    open(OUT, "> $resfile") or die "cannot write $resfile\n";
#    print OUT join("\t", @head, "#R", "Replicate"), "\n";
    print OUT join("\t", $head[0], "#R", "Replicate"), "\n";
    for (my $j=0; $j<@res; $j++) {
	my $name = $tbl->[$res[$j]]->[0];
	my $repl = $cnt{$name};
	my $nr = 0;
	$nr++ while($name =~ m/R/g);
#	print OUT join("\t", @{$tbl->[$res[$j]]}, $nr, $repl), "\n";
	print OUT join("\t", $tbl->[$res[$j]]->[0], $nr, $repl), "\n";
    }
    close(OUT);
}

