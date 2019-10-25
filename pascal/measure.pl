#!/usr/bin/perl
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::General qw(FileTester GetBondLength GetAngle GetTorsion);

sub init;
sub getAtoms;

my ($bgfFile, $selection, $geomType);
my ($atmList, $ATOMS, $geomVal, $getGeomVal);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, undef) = GetBGFFileInfo($bgfFile,0);
$atmList = getAtoms($ATOMS, $selection);
print "Done\n";
$geomVal = $getGeomVal->(@{ $atmList });
print "$geomType = $geomVal\n";

sub getAtoms {
    my ($atoms, $atmSelect)  = @_;
    my ($i, @ATMLIST, $selectName);

    for $i (@{ $atmSelect }) {
	die "ERROR: Atom $i is not a valid atom!\n" if (! exists($atoms->{$i}));
	push @ATMLIST, $atoms->{$i};
	$selectName .= "$i ($atoms->{$i}{FFTYPE}) ";
    }
    $geomType = $selectName . $geomType;
    return \@ATMLIST;
}    

sub init {
    my (%OPTS, $atm);
    getopt('bta',\%OPTS);
    for ("b", "a") {
	die "usage: $0 -b bgf file -a \"atom numbers\" -t [measurement type (bond|angle|torsion)]\n"
	    if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($bgfFile, $atm, $geomType) = ($OPTS{b}, $OPTS{a}, $OPTS{t});
    FileTester($bgfFile);
    while ($atm =~ /(\d+)/g) {
	push @{ $selection }, $1;
    }
    die "ERROR: Expected integer for atoms! Got \"$atmList\"\n" if (! $selection);
    die "ERROR: Expected more than 1 valid atom number!\n" if ($#{ $selection } == 0);
    if ($#{ $selection } > 3) {
	print "warning: max number of atoms should be 4..taking first 4...";
	while ($#{ $selection } > 3) {
	    pop @{ $selection };
	}
    }
    
    if ($#{ $selection } == 1) { #bond
	print "will compute bond..." if (defined($geomType) and $geomType !~ /^bond$/i);
	$geomType = "Bond";
	$getGeomVal = \&GetBondLength;
    } elsif ($#{ $selection } == 2) { 
	print "will compute angle..." if (defined($geomType) and $geomType !~ /^angle$/i);
	$geomType = "Angle";
	$getGeomVal = \&GetAngle;
    } else {
	print "will compute torsion..." if (defined($geomType) and $geomType !~ /^(torsion|dihedral)$/i);
	$geomType = "Torsion";
	$getGeomVal = \&GetTorsion;
    }
    print "Done\n";
}
