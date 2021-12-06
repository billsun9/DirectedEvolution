#!/usr/bin/perl

use strict;
use warnings;

# usage:
# cd /home/saito/work/160215_aiseq/data/Umetsu_GFP
# ./this_program | sort | uniq > Umetsu_ABLys.csv

my $OriDir = "original";
my $InFile = "$OriDir/Umetsu_GFP.csv";

# GFP = "SSHT"
# YFP = "GAYF"
my $tmpseq;
my @filelist;

# 0,1,2,3 = 65,72,77,203
my %p2i;
$p2i{"65"} =  0;
$p2i{"72"} =  1;
$p2i{"77"} =  2;
$p2i{"203"} = 3;

# $uniqseq{"SSHT"} = 1;
my %uniqseq;


print join(",", ("Sequence","Intensity","Change")), "\n";

$tmpseq = "SSHT";
@filelist = ($InFile);
for (my $i=0; $i<@filelist; $i++) {
    my $file = $filelist[$i];
    print STDERR "processing $file\n";
    open(IN, $file) or die "cannot read $file\n";
    while (my $line=<IN>) {
	chomp $line;
	$line =~ s/\,+$//;
	my $seq = $tmpseq;
	my ($aa, $ins, $chg) = split(/\,/, $line);
	$aa = uc($aa);
	$aa =~ s/\s+//g;
	if ($aa eq "GFPHIS") {
	    $seq = "SSHT";
	}
	elsif ($aa eq "YFPHIS") {
	    $seq = "GAYF";
	}
	elsif ($aa =~ /^(\S)(\d+)(\S)$/) {
	    my ($a, $p, $b) = ($1, $2, $3);
	    exists($p2i{$p}) or die "wrong mutation position $aa\n";
	    substr($seq, $p2i{$p}, 1) eq $a or die "wrong mutation character $aa\n";
	    substr($seq, $p2i{$p}, 1) = $b;
	}
	elsif ($aa =~ /^(\S{4})$/) {
	    $seq = $aa;
	}
	else {
	    die "never reach here $aa\n";
	}

	if (exists($uniqseq{$seq})) { 
# So far, all multiple entries are intensity=0, change=na. 
# Therefore, we simply skip them. 
	    print STDERR "warning: multiple entries $seq\n";
	    next;
	}

	&post_process(\$ins, \$chg, $seq) and next;
	print join(",", ($seq, $ins, $chg)), "\n";
	$uniqseq{$seq} = 1; 
    }
    close(IN);
}


# The measurements of intensity and change are not independent. 
# When intensity is very weak, the measurement of change is not reliable. 
# In this dataset, very weak intensity is regarded as 0, 
# and the corresponding change is regarded as a missing value (na). 
# We retain these entries as nagative data (0, 0) for machine learning. 
# In addition, we normalize intensity and change values 
# so that GFP has intensity=1 and change=1.
sub post_process {
    my ($ins_ref, $chg_ref, $seq) = @_;
    my ($ins_gfp, $chg_gfp) = (3065.668203, 0.750478377);

    $$ins_ref = $$ins_ref / $ins_gfp;
    $$chg_ref eq "na" and $$chg_ref=0;
    $$chg_ref = $$chg_ref / $chg_gfp;

    return 0;
}
