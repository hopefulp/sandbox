#!/usr/bin/perl
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::AMBER qw(parseAmberLib);
use Packages::FileFormats qw(GetBGFFileInfo createMOL2);
use Packages::General qw(FileTester);

sub init;
sub typeAtoms;
sub numerically { ($a<=>$b); }

my ($bgfFile, $mol2File, $LIBS);
my ($ATOMS, $BONDS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, undef) = GetBGFFileInfo($bgfFile);
print "Done\nTyping atoms...";
&typeAtoms($ATOMS, $BONDS, $LIBS);
print "Done\nCreating MOL2 file $mol2File...";
&createMOL2($ATOMS, $BONDS, $mol2File);
print "Done\n";

sub typeAtoms {
    my ($atoms, $bonds, $libs) = @_;
    my ($i, $atmName, $resName);

    for $i (sort numerically keys %{ $atoms }) {
	($atmName, $resName) = ($atoms->{$i}{ATMNAME}, uc($atoms->{$i}{RESNAME}));
	if (exists($libs->{$resName})) {
	    if (exists($libs->{$resName}{atoms}{$atmName})) {
		$atoms->{$i}{TYPED} = 1;
	    }
	}
    }


}

sub init {
    my (%OPTS, $libStr, @tmp, $parms, $i);
    getopt('bml',\%OPTS);
    die "usage: $0 -b bgf file -l (\"amber lib1 amber lib2...\") -m (mol2 name)\n"
	if (! exists($OPTS{b}));
    print "Initializing...";
    ($bgfFile, $libStr, $mol2File) = ($OPTS{b}, $OPTS{l}, $OPTS{m});
    FileTester($bgfFile);
    $libStr = "/exec/amber9/dat/leap/lib/all_amino02.lib /exec/amber9/dat/leap/lib/all_nucleic02.lib"
	if (! defined($libStr));
    while ($libStr =~ /(\S+)/g) {
	push @tmp, $1 if (-e $1 and -r $1 and -T $1);
    }
    die "ERROR: No valid amber libs found while search $libStr!\n" if (! @tmp);
    if (! defined($mol2File)) {
	$mol2File = basename($bgfFile);
	$mol2File =~ s/\.\w+$//;
	$mol2File .= "_amber.mol2";
    }
    print "Loading amber libs...";
    for $i (@tmp) {
	$LIBS = &parseAmberLib($i, $LIBS);
    }
    print "Done\n";
}
