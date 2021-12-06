#!/usr/bin/perl

use strict;
use warnings;

# usage:
# ./this_program Umetsu_GFP.csv psf seqlen > Umetsu_GFP_psf.csv
# ./this_program Umetsu_GFP.csv VHSE psf seqlen > Umetsu_GFP_VHSE_psf.csv

our $RS = "/home/saito/work/160215_aiseq/resource/aafeat";
our $RS2 = "/home/saito/work/160215_aiseq/data/Umetsu_GFP/psf";

# from [Westen et al, J Cheminform, 2013]
my %FeatFile = (
# BLOSUM.txt is extracted from Table 5 in [Georgiev, J Comput Biol, 2009], 
# which corresponds to features derived by eigenvalue decomposition of BLOSUM. 
# In [Westen et al, J Cheminform, 2013], they describe  
# > a descriptor set based on a VARIMAX analysis of physicochemical properties 
# > which were subsequently converted to indices based on the BLOSUM62
# which implies they use Table 7 in [Georgiev, J Comput Biol, 2009]. 
# However, by inspecting the supplementary file in [Westen et al, J Cheminform, 2013], 
# I realized that they actually used Table 5 rather than Table 7. 
"BLOSUM" => "$RS/BLOSUM.txt", 
"FASGAI" => "$RS/FASGAI.txt",
"MS-WHIM" => "$RS/MS-WHIM.txt",
"T-scale" => "$RS/T-scale.txt",
"ST-scale" => "$RS/ST-scale.txt",
"Z-scale" => "$RS/Z-scale.txt",
"VHSE" => "$RS/VHSE.txt",
# Features proposed in [Westen et al, J Cheminform, 2013]
"ProtFP" => "$RS/ProtFP.txt",
"ProtFP-Feature" => "$RS/ProtFP-Feature.txt",
);
# Position specific feature
my %PsfFile = (
# currently no psf
);

# $FeatValue->{"A"} = [0.01, -0.24, 0.35];
my $FeatValue = {};
# $PsfValue->[0]->{"A"} = [0.01, -0.24, 0.35];
my $PsfValue = [];

my $InFile = shift(@ARGV);
my $SeqLen = pop(@ARGV);
my @FeatList = @ARGV;

&init_feat($FeatValue);
for (my $i=0; $i<@FeatList; $i++) {
    exists($FeatFile{$FeatList[$i]}) or next;
    &load_feat($FeatValue, $FeatFile{$FeatList[$i]});
    print STDERR "load feature $FeatFile{$FeatList[$i]}\n";
}
foreach my $aa (sort keys %$FeatValue) {
    my @feat = @{$FeatValue->{$aa}};
    print STDERR "$aa feature ", join(",", @feat), "\n";
}
&init_psf($PsfValue,$SeqLen);
for (my $i=0; $i<@FeatList; $i++) {
    exists($PsfFile{$FeatList[$i]}) or next;
    &load_psf($PsfValue, $SeqLen, $PsfFile{$FeatList[$i]});
    print STDERR "load psf $PsfFile{$FeatList[$i]}\n";
}
for (my $i=0; $i<$SeqLen; $i++) {
    foreach my $aa (sort keys %{$PsfValue->[$i]}) {
	my @psf = @{$PsfValue->[$i]->{$aa}};
	print STDERR "position $i $aa psf ", join(",", @psf), "\n";
    }
}


my $tbl;
my %hidx;

open(IN, $InFile) or die "cannot read $InFile\n";
while (my $line=<IN>) {
    chomp $line;
    my @cells = split(/\,/, $line);
    for (my $i=0; $i<@cells; $i++) {
	$cells[$i] =~ /^\"(.*)\"$/ and $cells[$i] = $1;
    }

    if (!%hidx) {
	for (my $i=0; $i<@cells; $i++) {
	    $hidx{$cells[$i]} = $i;
	    print STDERR "load header entry $cells[$i]\n";
	}
	next;
    }

    my $name; 
    my @feat;
    my $pred;
    $name = $cells[$hidx{"Sequence"}];
    push(@feat, &get_feat($cells[$hidx{"Sequence"}], $FeatValue, $PsfValue, $SeqLen));
    $pred = &SigmoidProduct_pred($cells[$hidx{"Intensity"}], $cells[$hidx{"Change"}]);

    push(@$tbl, [$name, @feat, $pred]);
}
close(IN);

# Since combo accepts negative socres only, 
# we need to transform raw scores into negative values. 
# In this case, raw scores are normalized into the range 0-1 
# by using sigmoldal functions.
# Therefore, we do the following transformation 
# transformed_score = raw_score - 1.
my $max = 1;
my $idx = @{$tbl->[0]}-1;
for (my $i=0; $i<@$tbl; $i++) {
    $tbl->[$i]->[$idx] -= $max;
}
# Alternatively, we might employ the following transformations: 
# transformed_score = raw_score - max(observed_score)
# transformed_score = -exp(-raw_score)

my @head;
push(@head, "Sequence");
for (my $i=0; $i<@{$tbl->[0]}-2; $i++) {
    push(@head, "f");
}
push(@head, "Score");

print join(",", @head), "\n";
for (my $i=0; $i<@$tbl; $i++) {
    print join(",", @{$tbl->[$i]}), "\n";
}


sub init_feat {
    my ($FeatValue) = @_;
    foreach my $aa ("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V") {
	$FeatValue->{$aa} = [];
    }    
}
sub init_psf {
    my ($PsfValue, $SeqLen) = @_;
    for (my $i=0; $i<$SeqLen; $i++) {
	$PsfValue->[$i] = {};
	foreach my $aa ("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V") {
	    $PsfValue->[$i]->{$aa} = [];
	}
    }
}

sub load_feat {
    my ($FeatValue, $file) = @_;

    open(IN, $file) or die "cannot read $file\n";
    while (my $line=<IN>) {
	chomp $line;
	my @cells = split(/\t/, $line);
	my $aa = shift(@cells);
	my @feat = @cells;
	push(@{$FeatValue->{$aa}}, @feat);
    }
    close(IN);
}
sub load_psf {
    my ($PsfValue, $SeqLen, $file) = @_;

    for (my $i=0; $i<$SeqLen; $i++) {
	open(IN, "$file.$i.txt") or die "cannot read $file.$i.txt\n";
	while (my $line=<IN>) {
	    chomp $line;
	    my @cells = split(/\t/, $line);
	    my $aa = shift(@cells);
	    my @psf = @cells;
	    push(@{$PsfValue->[$i]->{$aa}}, @psf);
	}
	close(IN);
    }
}


sub get_feat {
    my ($sequence, $FeatValue, $PsfValue, $SeqLen) = @_;
    my @feat;

    length($sequence)==$SeqLen or die "wrong sequence length $sequence $SeqLen\n";

    for (my $i=0; $i<length($sequence); $i++) {
	my $aa = substr($sequence, $i, 1);
	push(@feat, @{$FeatValue->{$aa}});
	push(@feat, @{$PsfValue->[$i]->{$aa}});
    }

    return @feat;
}

sub sigmoid {
    my ($x) = @_;
    return 1.0 / (1.0 + exp(-$x));
}

sub SigmoidProduct_pred {
    my ($ins, $chg) = @_;
    my $thsh = 1.0;
    my $pred = &sigmoid($ins-$thsh) * &sigmoid($chg-$thsh);

    print STDERR "$ins $chg pred $pred\n";
    return $pred;
}
