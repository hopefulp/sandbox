#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetPDBFileInfo createPDB);
use Packages::General qw(FileTester);

sub init;
sub removeHydrogens;

my ($pdbFile, $saveFile);
my ($ATOMS, $BONDS);

$|++;
&init;
print "Parsing PDBFile $pdbFile...";
($ATOMS, $BONDS) = GetPDBFileInfo($pdbFile);
print "Done\nRemoving hydrogens...";
&removeHydrogens($ATOMS);
print "Done\nSaving PDB file $saveFile...";
createPDB($ATOMS, undef, $saveFile);
print "Done\n";

sub removeHydrogens {
    my ($atoms) = $_[0];
    my ($i);

    for $i (keys %{ $atoms }) {
	delete($atoms->{$i}) if ($atoms->{$i}{ATMNAME} =~ /^\s*\d*H/);
    }
}
    
sub init {
    my (%OPTS);
    getopt('ps',\%OPTS);
    die "usage: $0 -p pdb file -s [save name]\n" if (! exists($OPTS{p}));
    print "Initializing...";
    ($pdbFile, $saveFile) = ($OPTS{p}, $OPTS{s});
    FileTester($pdbFile);
    if (! defined($saveFile)) {
	$saveFile = basename($pdbFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "nopdb.pdb";
    }
    print "Done\n";
}
