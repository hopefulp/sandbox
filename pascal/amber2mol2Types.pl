#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetMOL2FileInfo createMOL2 AMBER2MOL2Types);
use Packages::General qw(FileTester);

sub init;

my ($mol2File, $saveFile);

$|++;
&init;
my ($ATOMS, $BONDS) = GetMOL2FileInfo($mol2File, 0);
print "Done\nConverting atom types...";
delete $ATOMS->{HEADER};
&AMBER2MOL2Types($ATOMS);
print "Done\nCreating MOL2 file $saveFile...";
createMOL2($ATOMS, $BONDS, $saveFile, 1);
print "Done\n";

sub init {
    my (%OPTS);
    getopt('ms',\%OPTS);
    die "usage: $0 -m mol2file -s [save file]\n" if (! exists($OPTS{m}));
    print "Initializing...";
    ($mol2File, $saveFile) = ($OPTS{m}, $OPTS{s});
    FileTester($mol2File);
    if (! defined($saveFile)) {
	$saveFile = basename($mol2File);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_converted.mol2";
    }
    print "Done\n";
}
	
