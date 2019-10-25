#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::General qw(FileTester);
use Packages::FileFormats qw(createBGF addHeader GetBGFFileInfo);

sub init;
sub updateAtomTypes;
sub updateResData;
sub updateFieldData;

my ($resFile, $saveFile, $bgfFile, $FIELDS);
my ($ATOMS, $BONDS, $HEADERS, $resATOMS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nParsing residue library file $resFile...";
($resATOMS, undef, undef) = GetBGFFileInfo($resFile);
$resATOMS = updateResData($resATOMS);
&updateFieldData($ATOMS, $FIELDS);
print "Done\nUpdating atom info based on library...";
&updateAtomTypes($ATOMS, $resATOMS, $FIELDS);
print "Done\nCreating BGF file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub updateFieldData {
    my ($atoms, $fields) = $_[0];
    my ($i, @tmp, $fatom, $fList);

    @tmp = keys %{ $atoms };
    $fatom = $atoms->{ $tmp[0] };

    for $i (keys %{ $fields }) {
	$fList .= "$i ";
	delete $fields->{$i} if (! exists($fatom->{$i}));
    }

    die "ERROR: field \"$fList\" is invalid!\n" if (! scalar(keys %{ $fields }) == 0);
}

sub updateAtomTypes {
    my ($atoms, $resInfo, $fields) = @_;
    my ($i, $j, $resName, $atmName);

    for $i (keys %{ $atoms }) {
	($resName, $atmName) = ($atoms->{$i}{RESNAME}, $atoms->{$i}{ATMNAME});
        $resName =~ s/\s//g;
        $atmName =~ s/\s//g;	
	if (exists($resInfo->{$resName}) and exists($resInfo->{$resName}{$atmName})) {
	    for $j (keys %{ $fields }) {
		$atoms->{$i}{$j} = $resInfo->{$resName}{$atmName}{$j};
	    }
	}
    }
}

sub updateResData {
    my ($resInfo) = $_[0];
    my (%RES, $i, $resName, $atmName);

    for $i (keys %{ $resInfo }) {
	$resName = $resInfo->{$i}{RESNAME};
	$atmName = $resInfo->{$i}{ATMNAME};
	$resName =~ s/\s//g;
	$atmName =~ s/\s//g;
	%{ $RES{$resName}{$atmName} } = %{ $resInfo->{$i} };
    }
    
    return \%RES;
}

sub init {
    my (%OPTS, $fList);
    getopt('brsf',\%OPTS);
    
    for ("b", "r") {
	die "usage: $0 -b bgf structure file -r residue library (bgf) -s (save name (optional) -f (fields = fftype)\n"
	    if (! exists($OPTS{$_}));
    }

    ($bgfFile, $saveFile, $resFile, $fList) = ($OPTS{b}, $OPTS{s}, $OPTS{r}, $OPTS{f});
    print "Initializing...";
    FileTester($bgfFile);
    FileTester($resFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_typed.bgf";
    }

    $fList = "FFTYPE" if (! defined($fList));
    while ($fList =~ /(\S+)/g) {
	$FIELDS->{uc($1)} = 1;
    }
    print "Done\n";
}
