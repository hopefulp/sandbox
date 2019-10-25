#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use warnings;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use Packages::General qw(FileTester);
use Packages::ManipAtoms qw(GetAtmList GetMols);


sub init;
sub parseAtomPosFile;
sub saveSlabDataAsChains;
sub numerically { ($a<=>$b); }

my ($bgfFile, $atmPosFile, $saveFile, $num, $shellDist);
my ($ATOMS, $BONDS, $HEADERS, $MOLS);

$|++;
&init;
print "Parsing bgf file...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$MOLS = GetMols($ATOMS, $BONDS);
print "Done\nParsing atom position file $atmPosFile...";
&parseAtomPosFile($ATOMS, $atmPosFile);
print "Done\nCreating $num groups...";
&saveSlabDataAsChains($ATOMS, $MOLS, $num, $shellDist);
print "Done\nWriting BGF file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub saveSlabDataAsChains {
    my ($atoms, $mols, $tot, $dists) = @_;
    my ($dens, $i, @list, $min, $max, $index);
    my ($j, $count, $inc, $k, $avg);

    for $i (keys %{ $atoms }) {
	$dens->{$atoms->{$i}{POS}{AVG}}++;
	$min = $atoms->{$i}{POS}{AVG} if (! defined($min) or $atoms->{$i}{POS}{AVG} < $min);
	$max = $atoms->{$i}{POS}{AVG} if (! defined($max) or $atoms->{$i}{POS}{AVG} > $max);
    }
    @list = sort numerically keys %{ $dens };
    if(defined($dists) and $list[$#list] < $dists->[$#{ $dists }]) {
	$list[$#list+1] = $dists->[$#{ $dists }]; #if the user specified more distances than we find differences
						    #make the greatest entry in the list the largest specified distance
    }

    if (! defined($dists)) {
	#need to make equal intervals
	$num = $#list+1 if(($#list+1)<$num);
	$inc = (($max-$min)/$num); # increment based on distance distributions
	$dists->[0] = $min + 0.9*$inc;
	$dists->[$num-1] = $max;
	for $i (1 .. ($num-2)) {
	    $dists->[$i] = $min+$inc*($i+0.9);
	}
    }

    for $i (keys %{ $mols }) {
	$avg = $count = 0;
	for $j (keys %{ $mols->{$i}{MEMBERS} }) {
	    $count++;
	    $avg += $atoms->{$j}{POS}{AVG};
	}
	$avg = sprintf("%.0f",$avg/$count);
        $index = 64+$num;
	for $k (0 .. ($num-1)) {
	    $index = 64 + $k+1;
	    last if ($avg < $dists->[$k]);
	}
	for $j (keys %{ $mols->{$i}{MEMBERS} }) {
	    $atoms->{$j}{CHAIN} = chr($index);
	}
    }
}

sub parseAtomPosFile {
    my ($atoms, $infile) = @_;
    my ($count, $tot);

    $tot = scalar(keys %{ $ATOMS });
    open ATOMPOSFILE, $infile or die "ERROR: Cannot open $infile: $!\n";
    while(<ATOMPOSFILE>) {
	chomp;
	if ($_ =~ /^\s+(\d+)\s+(\d+)\s+\-?(\d+\.\d+)/) {
	    next if (! exists($atoms->{$1}));
	    $atoms->{$1}{POS}{AVG} = $2;
	    $atoms->{$1}{POS}{STDEV} = $3;
	    $count++;
	}
    }
    close ATOMPOSFILE;
    die "ERROR: $infile has $count entries while bgf file has $tot atoms!\n"
	if($count != $tot);
}

sub init {
    my (%OPTS, $select, @tmp, $shellDistStr);

    getopt('basnd',\%OPTS);
    die "usage: $0 -b bgf file -a atom postion file -n [numgroups = 10] -d [shell distances] -s [save name]\n"
	if(! defined($OPTS{b}) or ! defined($OPTS{a}));

    print "Initializing...";
    ($bgfFile, $atmPosFile, $saveFile, $num, $shellDistStr) = ($OPTS{b}, $OPTS{a}, $OPTS{s}, $OPTS{n}, $OPTS{d});
    FileTester($bgfFile);
    FileTester($atmPosFile);
    $num = 10;
    die "ERROR: Expected positive integer for $num, got \"$num\"\n" if ($num !~ /^\d+$/ or !$num);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+//;
	$saveFile .= "_slab.bgf";
    }
    undef $shellDistStr if(defined($shellDistStr) and $shellDistStr !~ /^\d+\.?\d*/);
    if (defined($shellDistStr)) {
	$num = 0;
	while ($shellDistStr =~ /(\d+\.?\d*)/g) {
	    push @{ $shellDist }, $1;
	    $num++;
	}
    }
    print "Done\n";
}
