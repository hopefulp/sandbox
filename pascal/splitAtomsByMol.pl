#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use warnings;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use Packages::General qw(FileTester);
use Packages::ManipAtoms qw(SelectAtoms BuildAtomSelectionString GetMols);

sub init;
sub splitAtoms;
sub showUsage;

my ($bgfFile, $saveFile, $SELECT, $atomSelect, $field);
my ($ATOMS, $BONDS, $HEADERS, $MOLS);

$|++;
&init;
print "Parsing bgf file...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$SELECT = SelectAtoms($atomSelect, $ATOMS);
die "No selected atoms found in bgf file\n" if (! keys %{ $SELECT });
print "Done\nDetermining molecules based on connectivity...";
$MOLS = &GetMols($ATOMS,$BONDS,$SELECT);
print "Done\nSpliting atoms by molecules...";
&splitAtoms($ATOMS, $MOLS, $SELECT);
print "Done\nCreating bgf file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub splitAtoms {
    my ($atoms, $mols, $selection) = @_;
    my ($i);

    for $i (keys %{ $selection }) {
	if ($field eq "CHAIN") {
	    $atoms->{$i}{$field} = chr((64 + ${ $atoms->{$i}{MOLECULEID} }));
	    $atoms->{$i}{$field} = "X" if (${ $atoms->{$i}{MOLECULEID} } > 4);
	} else {
	    $atoms->{$i}{$field} = ${ $atoms->{$i}{MOLECULEID} };
	}
    }
}

sub init {
    my (%OPTS, $select, $usage);

    getopt('bafs',\%OPTS);
    $usage = showUsage();
    die "$usage\n" if (! defined($OPTS{b}));
    print "Initializing...";
    ($bgfFile, $saveFile, $select, $field) = ($OPTS{b}, $OPTS{s}, $OPTS{a}, $OPTS{f});
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+//;
	$saveFile .= "_mod.bgf";
    }
    
    $select = "index > 0" if (! defined($select));
    $atomSelect = BuildAtomSelectionString($select);
    $field = "CHAIN" if (! defined($field) or $field !~ /^CHAIN|RESNUM|RESID|RES$/i);
    $field = uc($field);
    $field = "RESNUM" if ($field =~ /RES/);

    print "Done\n";
}

sub showUsage {
    my ($usage) = <<DATA;
usage: $0 -b bgf file -a (atom selection) -f (field) -s (save name)
Options:
	-b bgf file (Required): location of bgf structure file
	-a atom selection (optional): any valid field (fftype, charge, resname etc) expression
		e.g. -a "resname ne 'WAT' and resnum > 10". Default: all atoms
	-f field (optional): field to use to split molecules. Either chain or resid (default)
	-s save name: name of new bgf file
DATA
    return $usage;
}

