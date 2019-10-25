#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use Packages::General qw(FileTester);
use Packages::FileFormats qw(GetBGFFileInfo createMOL2);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

my ($ATOMS, $BONDS, $HEADERS, $bgfFile, $mol2File, $verbose);
$|++;
&init;
print "Parsing BGF file $bgfFile..." if ($verbose);
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1, $verbose);
print "Done\nCreating MOL2 file $mol2File..." if ($verbose);
&createMOL2($ATOMS, $BONDS, $mol2File, $forAmber);
print "Done\n" if ($verbose);

sub init {
    my (%OPTS);
    getopt('bamv',\%OPTS);
    ($bgfFile, $mol2File, $forAmber, $verbose) = ($OPTS{b}, $OPTS{m}, $OPTS{a}, $OPTS{v});
    die "usage: $0 -b bgfFile -m [mol2file] -a [mol24amber?=no] -v [verbose=yes]\n" if (! defined($bgfFile));
    $verbose = 1 if (! defined($verbose) or $verbose !~ /^0|no$/i);
    $verbose = 0 if ($verbose =~ /^0|no$/i);
    print "Initializing..." if ($verbose);
    FileTester($bgfFile);
    if (! defined($mol2File)) {
	$mol2File = basename($bgfFile);
	$mol2File =~ s/\.\w+$/\.mol2/;
    }
    $forAmber = 1 if (! defined($forAmber));
    if ($forAmber =~ /1|yes/i) {
	$forAmber = 1;
    } else {
	$forAmber = 0;
    }
    print "Done\n" if ($verbose);
}
