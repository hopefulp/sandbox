#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo createBGF addHeader);
use Packages::LAMMPS qw (ReadDataFile);
use Packages::General qw (FileTester);
use File::Basename;
use Getopt::Std qw(getopt);

sub init;
sub updateCoordinates;
sub updateBonds;

my ($dataFile, $bgfFile, $uType, $saveFile);
my ($ATOMS, $BONDS, $HEADERS, $DATA);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nParsing LAMMPS data file $dataFile...";
$DATA = ReadDataFile($dataFile);
print "Done\nCreating updated BGF file $saveFile...";
updateCoordinates($ATOMS, $DATA->{ATOMS}) if ($uType == 1);
updateBonds($BONDS, $DATA->{BONDS}) if ($uType == 2);
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub init {
    my (%OPTS);

    getopt('bdsu',\%OPTS);

    for ("b", "d") {
	die "usage: $0 -b bgf file -d datafile -o (update: 1=atom data, 2=bond data) -s (savename)\n"
	    if (! exists($OPTS{$_}));
    }
 
    ($bgfFile, $dataFile, $uType, $saveFile) = ($OPTS{b}, $OPTS{d}, $OPTS{u}, $OPTS{s});
    print "Initializing...";   
    FileTester($bgfFile);
    FileTester($dataFile);
    $uType = 1;
    $uType = $1 if ($uType =~ /(1|2)/);
    if (! $saveFile) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//g;
	$saveFile .= "_updated.bgf";
    }
    print "Done\n";
}

sub updateCoordinates {
    my ($ATOMS, $COORDS) = @_;
    my ($i, $dim, $j, @tmp, @tmp1);

    @tmp =  keys %{ $ATOMS };
    @tmp1 = keys %{ $COORDS };

    die "ERROR: Inconsistent atom count in BGF and DATA files... aborting\n"
	if ($#tmp != $#tmp1);

    for $i (@tmp) {
	$j = 3;
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    $ATOMS->{$i}{$dim} = $COORDS->{$i}{$j};
	    $j++;
	}
    }
}

sub updateBonds {
    my ($BONDS, $CONS) = @_;
    my ($i, $j);

    $BONDS = ();

    for $i (sort {$a<=>$b} keys %{ $CONS }) {
	for $j (keys %{ $CONS }) {
	    $CONS->{$j}{$i} = 1;
	    push @{ $BONDS->{$i} }, $j;
	}
    }
}
