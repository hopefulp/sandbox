#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use Packages::General qw(GetSelections FileTester);
use Packages::ManipAtoms qw(GetAtmList);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub updateResNum;
sub numerically { ($a<=>$b); }

my ($saveFile, $bgfFile, $selection);
my ($ATOMS, $BONDS, $HEADERS, $SELECT, $SELECTATMS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nSelecting atoms...";
$SELECT = GetSelections($selection, 0);
$SELECTATMS = GetAtmList($SELECT, $ATOMS);
print "Done\nUpdating residue numbers to match molecules...";
&updateResNum($ATOMS, $SELECTATMS);
print "Done\nCreating BGF file $saveFile...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub updateResNum {
    my ($atoms, $selectAtoms) = @_;
    my ($i, $j, $resNum, $firstAtm);

    for $i (sort numerically keys %{ $selectAtoms }) {
	next if ($atoms->{$i}{UPDATED});
	if (! defined($resNum)) {
	    $firstAtm = $i;
	    if ($i == 1) {
		$resNum = 0;
	    } else {
		$resNum = $atoms->{$i -1}{RESNUM};
	    }
	}
	$resNum++;
	$atoms->{$i}{RESNUM} = $resNum;
	$atoms->{$i}{UPDATED} = 1;
	next if (! exists($atoms->{$i}{MOLECULE}));
	for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
	    $atoms->{$j}{RESNUM} = $resNum;
	    $atoms->{$j}{UPDATED} = 1;
	}
    }

    for $i  (keys %{ $atoms }) {
	next if ($atoms->{$i}{UPDATED} or $i <= $firstAtm);
	$atoms->{$i}{RESNUM} += $resNum;
    }
}
	    
sub init {
    my ($atomSelect, %OPTS);

    getopt('bsa',\%OPTS);
    for ("b", "a") {
	die "usage: $0 -b bgf file -a atom selection -s [save name]\n" if (! defined($OPTS{$_}));
    }

    print "Initializing...";
    ($bgfFile, $saveFile, $atomSelect) = ($OPTS{b}, $OPTS{s}, $OPTS{a});
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_mod.bgf";
    }

    if ($atomSelect =~ /\s+/) {
        @{ $selection } = split /\s+/, $atomSelect;
    } else {
        $selection->[0] = $atomSelect;
    }

    print "Done\n";
}
