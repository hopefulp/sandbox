#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use Packages::FileFormats qw(GetBGFFileInfo GetMOL2FileInfo GetPDBFileInfo addHeader
			     createBGF createPDB createMOL2 insertHeaderRemark);
use Packages::General qw(FileTester);
use Packages::ManipAtoms qw(SelectAtoms BuildAtomSelectionString GetMols);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub showUsage;
sub init;
sub updateAtomFields;
sub parseFieldList;
sub includeMolAtoms;

my ($FILES, $FIELDS, $selection, $SELECT, $saveType, $fieldStr, $molOpt);
my ($ATOMS, $BONDS, $HEADERS, $SELECTATMS, $saveFunc, $readFunc);

$|++;
&init;
print "Parsing Structure file $FILES->{STRUCTURE}...";
($ATOMS, $BONDS, $HEADERS) = $readFunc->($FILES->{STRUCTURE},1);
&GetMols($ATOMS, $BONDS);
if ($saveType =~ /mol2/) {
    $HEADERS = $ATOMS->{HEADER};
    delete $ATOMS->{HEADER};
}
&parseFieldList($ATOMS, $FIELDS);
print "Done\nParsing atom/residue selection...";
$SELECTATMS = SelectAtoms($selection, $ATOMS);
die "ERROR: Not a valid atom selection!\n" if (! keys %{ $SELECTATMS });
&includeMolAtoms($ATOMS, $SELECTATMS) if ($molOpt);
print "Done\nUpdating fields..";
&updateAtomFields($ATOMS, $SELECTATMS, $FIELDS);
print "Done\nCreating $FILES->{SAVE}...";
&insertHeaderRemark($HEADERS, "REMARK $FILES->{STRUCTURE} modified $fieldStr") 
	if ($saveType eq "bgf");
&addHeader($ATOMS, $HEADERS) if ($saveType eq "bgf");
$saveFunc->($ATOMS, $BONDS, $FILES->{SAVE});
print "Done\n";

sub includeMolAtoms {
    my ($atoms, $select) = @_;
    my ($i, $j);

    for $i (keys %{ $atoms }) {
	next if (! exists($select->{$i}));
	for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
	    $select->{$j} = 1;
	}
    }
}

sub updateAtomFields {
    my ($atoms, $atomSelection, $fields) = @_;
    my ($i, $j, $tmp);

    for $i (keys %{ $atomSelection }) {
	for $j (keys %{ $fields }) {
	    if ($fields->{$j}{MOD} eq ".") {
		$atoms->{$i}{$j} .= $fields->{$j}{VAL};
	    } elsif ($fields->{$j}{MOD} eq "+") {
		$atoms->{$i}{$j} += $fields->{$j}{VAL};
	    } elsif ($fields->{$j}{MOD} eq "-") {
		$atoms->{$i}{$j} -= $fields->{$j}{VAL};
	    } else {
		$atoms->{$i}{$j} = $fields->{$j}{VAL};
	    }
	}
    }
}

sub parseFieldList {
    my ($atoms, $fields) = @_;
    my ($i, $currAtm, $fieldList, $j);

    for $i (keys %{ $atoms }) {
	$currAtm = \%{ $atoms->{$i} };
	for $j (keys %{ $fields }) {
	    delete $fields->{$j} if (! exists($currAtm->{$j}));
	    $fieldList .= "$j ";
	}
	last;
    }

    die "ERROR: No valid fields found while search for ${fieldList} in bgf file\n"
	if (! keys %{ $fields });
}

sub init {
    my (%OPTS, $usage, $atomSelect, $tmp, $rec);

    getopt('safwtm',\%OPTS);

    $usage = &showUsage;
    for ("s","a","f") {
	die "$usage\n" if (! defined($OPTS{$_}));
    }
    print "Initializing...";
    ($FILES->{STRUCTURE}, $atomSelect, $fieldStr, $FILES->{SAVE}, $saveType, $molOpt) = 
	($OPTS{s}, $OPTS{a}, $OPTS{f}, $OPTS{w}, $OPTS{t}, $OPTS{m});
    FileTester($FILES->{STRUCTURE});
    $atomSelect = "*" if (! defined($atomSelect));
    $selection = BuildAtomSelectionString($atomSelect);

    while ($fieldStr =~ /(\S+)/g) {
	$tmp = $1;
	if ($tmp =~ /(\w+):(\+|\-|\.)?(.*)/) {
	    if (defined($2)) {
		$rec = (
			{
			    "MOD" => $2,
			    "VAL" => $3,
			}
			);
	    } else {
		$rec = (
			{
			    "MOD" => "",
			    "VAL" => $3,
			}
			);
	    }
	    $FIELDS->{$1} = $rec;
	}
    }

    $saveType = "bgf" if (! defined($saveType));
    $saveType = lc($saveType);
    $saveFunc = \&createBGF;
    $readFunc = \&GetBGFFileInfo;
    if ($saveType =~ /mol2/) {
	$saveFunc = \&createMOL2;
	$readFunc = \&GetMOL2FileInfo;
    } elsif ($saveType =~ /pdb/) {
	$saveFunc = \&createPDB;
	$readFunc = \&GetPDBFileInfo;
    }
    if (! defined($FILES->{SAVE})) {
        $FILES->{SAVE} = $FILES->{STRUCTURE};
        $FILES->{SAVE} =~ s/\.\w+$//;
	$FILES->{SAVE} .= "_mod.${saveType}";
    }

    die "ERROR: Invalid string found while searching \"FIELD\".\n" .
	"Expected field:(+|-|.)new_val. Got \"${fieldStr}\"!\n" if (! defined($FIELDS));
    $molOpt = 0 if (! defined($molOpt) or $molOpt !~ /1|yes/i);
    $molOpt = 1 if ($molOpt =~ /1|yes/i);
    print "Done\n";
}

sub showUsage {
    return "usage: $0 -s structure file -a atom selection -f field -m [include mol] -w [save name]\n" .
	"Options:\n\t-s structure file: Location of structure file\n" .
	"\t-a atom selection: : [:][I|N|T][a|r]x[-y:z]\n" .
	"\t-f field: field to adjust. Expected field:[+|-|.]new_val. Enclose multiple in quotes\n" .
	"\t-m include mol: (Optional) Include all atoms of molecule if atom is selected. Default=no\n" .
	"\t-w [save name]: (Optional) Name to save the file as\n" .
	"\t-t [save type]: (Optional) File formats (bgf (default), mol2 or pdb) to save file as\n";
}
