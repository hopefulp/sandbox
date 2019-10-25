#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::FileFormats qw(GetMSIFileInfo addHeader createBGF);
use Packages::General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;

my ($msiFile, $bgfFile);
my ($ATOMS, $BONDS, $HEADERS);

$|++;
&init;
print "Parsing MSI file $msiFile...";
($ATOMS, $BONDS, $HEADERS) = GetMSIFileInfo($msiFile, 1);
print "Done\nCreating BGF file $bgfFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $bgfFile);
print "Done\n";

sub init {
    my (%OPTS);
    
    getopt('mb',\%OPTS);
    die "usage: $0 -m msi file -b (bgf file)\n" if (! exists($OPTS{m}));
    print "Initializing...";
    ($msiFile, $bgfFile) = ($OPTS{m}, $OPTS{b});
    FileTester($msiFile);
    if (! defined($bgfFile)) {
	$bgfFile = basename($msiFile);
	$bgfFile =~ s/\.\w+$//;
	$bgfFile .= ".bgf";
    }
    print "Done\n";
}
