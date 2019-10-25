#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createHeaders createBGF addBoxToHeader GetBGFAtoms insertHeaderRemark);;
use Packages::BOX qw(GetBox);
use Packages::ManipAtoms qw(BuildAtomSelectionString SelectAtoms GetMols);
use Packages::General qw(FileTester CoM);

sub init;
sub getNewBox;
sub trimCell;
sub makeAtomsMols;

my ($bgfFile, $saveFile, $newCell, $startOrigin, $isMol, $cStr);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $MOLS, $selection, $SELECT);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$BOX = GetBox($ATOMS, undef, undef);
if ($isMol) {
    $MOLS = &GetMols($ATOMS, $BONDS);
} else {
    $MOLS = makeAtomsMols($ATOMS);
}
&getNewBox($BOX, $newCell, $startOrigin);
print "Parsing atom/residue selection...";
$SELECT = SelectAtoms($selection, $ATOMS);
print "Done\nRemoving atoms outside new box...";
&trimCell($ATOMS, $BONDS, $MOLS, $BOX, $SELECT);
($ATOMS, $BONDS) = GetBGFAtoms($ATOMS, $ATOMS, $BONDS);
print "Done\nCreating BGF file $saveFile...";
&insertHeaderRemark($HEADERS, "REMARK $bgfFile trimed $cStr");
&addBoxToHeader($HEADERS, $BOX);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub makeAtomsMols {
    my ($atoms) = $_[0];
    my ($i, %atomsMols);

    for $i (keys %{ $atoms }) {
	$atomsMols{$i}{$i} = 1;
    }

    return \%atomsMols;
}

sub trimCell {
    my ($atoms, $bonds, $mols, $box, $select) = @_;
    my ($i, $molCenter, $j, $currMol, $isOutside);

    for $i (keys %{ $mols }) {
	$currMol = ();
	for $j (keys %{ $mols->{$i}{MEMBERS} }) {
	    %{ $currMol->{$j} } = %{ $atoms->{$j} };
	}
	$molCenter = CoM($currMol);
	$isOutside = 0;
	for $j (keys %{ $box }) {
	    if (($molCenter->{$j . "COORD"} > $box->{$j}{hi}) or
		($molCenter->{$j . "COORD"} < $box->{$j}{lo})) {
		$isOutside = 1;
		last;
	    }
	}
	if ($isOutside) {
	    for $j (keys %{ $mols->{$i}{MEMBERS} }) {
		last if (! exists($select->{$j}));
		delete $atoms->{$j};
		delete $bonds->{$j};
	    }
	}
    }
}

sub getNewBox {
    my ($oldBox, $newBox, $startPoint) = @_;
    my ($i, $offset, $isValid);

    $isValid = 0;
    for $i (keys %{ $oldBox }) {
	$offset = ($oldBox->{$i}{len} - $newBox->{$i})/2;
	if ($offset > 0) {
	    $isValid = 1;
	    $oldBox->{$i}{len} = $newBox->{$i};
	    if ($startPoint == 1) { # start at origin
		$oldBox->{$i}{hi} = $newBox->{$i};
		$oldBox->{$i}{lo} = 0;
	    } elsif ($startPoint == 2) { #start at min
		$oldBox->{$i}{hi} -= $offset*2;
	    } elsif ($startPoint == 3) { #start at max
		$oldBox->{$i}{lo} += $offset*2;
	    } else { # symmetric cut from box center
		$oldBox->{$i}{hi} -= $offset;
		$oldBox->{$i}{lo} += $offset;
	    }
	}
    }
    
    die "ERROR: New cell is larger than old cell!\n" if (! $isValid);

}

sub init {
    my (%OPTS, $atomStr, $usage);

    getopt('bcsoma',\%OPTS);
    $usage = showUsage();
    for ("b", "c") {
	die "$usage" if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($bgfFile, $cStr, $atomStr, $saveFile, $startOrigin, $isMol) = ($OPTS{b}, $OPTS{c}, $OPTS{a}, $OPTS{s}, $OPTS{o}, $OPTS{m});
    FileTester($bgfFile);
    $atomStr = "index > 0" if (! defined($atomStr));
    $selection = BuildAtomSelectionString($atomStr);

    $startOrigin = 0 if (! defined($startOrigin));
    $startOrigin = 1 if ($startOrigin =~ /^(1|yes)/i);
    $isMol = 1 if (! defined($isMol) or $isMol !~ /0|no/i);
    if ($cStr !~ /(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)/) {
        die "ERROR: Expected integers for x,y and z cell length. Got \"$cStr\"\n";
    } else {
        $newCell = (
		    {
			"X"      => $1,
			"Y"      => $2,
			"Z"      => $3,
			"STRING" => "${1}x${2}x${3}",
		    }
		    );
    }
                                                                                                                           
    if (! defined($saveFile)) {
        $saveFile = basename ($bgfFile);
        $saveFile =~ s/\.\w+$//;
        $saveFile .= "_trim_" . $newCell->{STRING} . ".bgf";
    }
    print "Done\n";
}

sub showUsage {
    my ($usage) = <<DATA;
usage: $0 -b bgf file -c \"x y z\" new cell length -a (atom selection) -s (save name) -o (origin) -m (keep together)
Options:
	-b bgf file (required): location of bgf structure file
	-c new cell length (required): expected "xx.xx yy.yy zz.zz" of new cell. must be larger than current cell
	-a atom selection (optional): atoms to consider when shrinking cell. Any valid field (fftype, charge etc)
		expression. e.g. -a "resname eq 'WAT' and resnum > 20". Default: all atoms
	-s save name (optional): name of new bgf file
	-m keep together (optional): flag to delete molecules instead of atoms. Default: 0
	-o origin (optional): 1: (0, 0, 0), 2: box min, 3: box max, 0 (default): box center
DATA
    return $usage;
}

