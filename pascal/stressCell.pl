#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createHeaders createBGF addBoxToHeader);
use Packages::BOX qw(GetBox);
use Packages::ManipAtoms qw(SplitAtomsByMol);
use Packages::General qw(FileTester CoM);

sub init;
sub getNewBox;
sub stressCell;

my ($bgfFile, $saveFile, $newCell);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $MOLS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, undef) = GetBGFFileInfo($bgfFile);
$BOX = GetBox($ATOMS, undef, undef);
$MOLS = SplitAtomsByMol($ATOMS, $BONDS);
&getNewBox($BOX, $newCell);
print "Done\nScaling cell by $newCell->{STRING}...";
&stressCell($ATOMS, $newCell);
print "Done\nCreating BGF file $saveFile...";
$HEADERS = createHeaders(undef, $saveFile);
&addBoxToHeader($HEADERS, $BOX);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub stressCell {
    my ($atoms, $cellScaling) = @_;
    my ($i, $j, $tot, $dimTotal, $dimScale, $cellExtrema);

    $tot = scalar keys %{ $atoms };
    for $i ("X", "Y", "Z") {
	$dimTotal = $cellExtrema = 0;
	for $j (keys %{ $atoms }) {
	    $dimTotal += $atoms->{$j}{$i . "COORD"};
	    $cellExtrema = $atoms->{$j}{$i . "COORD"} if (! $cellExtrema or $atoms->{$j}{$i . "COORD"} < $cellExtrema);
	}
	$dimScale = 1 - ($cellScaling->{$i}/(($dimTotal/$tot) - $cellExtrema));
	for $j (keys %{ $atoms }) {
	    $atoms->{$j}{$i . "COORD"} *= $dimScale;
	}
    }
}

sub getNewBox {
    my ($oldBox, $newBox) = @_;
    my ($i);

    for $i (keys %{ $oldBox }) {
	$oldBox->{$i}{hi} *= $newBox->{$i};
	$oldBox->{$i}{lo} *= $newBox->{$i};
	$oldBox->{$i}{len} *= $newBox->{$i};
    }
}

sub init {
    my (%OPTS, $cStr);
    getopt('bcs',\%OPTS);
    for ("b", "c") {
	die "usage: $0 -b bgf file -c \"x y z\" cell length scaling -s (save name)\n"
	    if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($bgfFile, $cStr, $saveFile) = ($OPTS{b}, $OPTS{c}, $OPTS{s});
    FileTester($bgfFile);
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
        $saveFile .= "_compress_" . $newCell->{STRING} . ".bgf";
    }
    print "Done\n";
}
