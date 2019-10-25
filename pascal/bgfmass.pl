#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::General qw(FileTester LoadFFs GetSelections);
use Packages::FileFormats qw(GetBGFFileInfo AddMass);
use Getopt::Std qw(getopt);
use Packages::ManipAtoms qw(GetAtmList);
use Packages::CERIUS2 qw(ReadFFs);

sub init;
sub calcMass;

my ($bgfFile, $ff);
my ($ATOMS, $SELECT);

$|++;
&init;
&calcMass($ATOMS);

sub init {
    my ($BONDS, $FFDATA, %OPTS, $FF, $select, $selection, $FFILES);

    getopt('bfa',\%OPTS);
    for ("b", "f") {
	die "usage: $0 -b bgf file -f force field -a (atom selection = all)\n"
	    if (! exists($OPTS{$_}));
    }

    print "Initializing...";
    ($bgfFile, $FF, $select) = ($OPTS{b}, $OPTS{f}, $OPTS{a});
    FileTester($bgfFile);
    $select =  "*" if (! defined($select));
    if ($select =~ /\s+/) {
        @{ $selection } = split /\s+/, $select;
    } else {
        $selection->[0] = $select;
    }

    print "Done\nParsing $bgfFile...";
    ($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
    print "Done\n";
    ($FFILES, undef) = ReadFFs($FF);
    $FFDATA = LoadFFs($FFILES);
    print "Adding Mass...";
    AddMass($ATOMS, $FFDATA);
    print "Done\nSelecting atoms...";
    $SELECT = GetSelections($selection, 0);
    $SELECT = GetAtmList($SELECT, $ATOMS);
    print "Done\n";
}

sub calcMass {
    my ($atoms) = $_[0];
    my ($i, $mass);

    $mass = 0;
    for $i (keys %{ $atoms }) {
	$mass += $atoms->{$i}{MASS};
    }
    
    print "Total Mass: $mass amu\n";
}
