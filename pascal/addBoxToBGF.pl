#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use Packages::General qw(FileTester);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use Packages::BOX qw(GetBox);
use Packages::CERIUS2 qw(parseCerius2FF);
use File::Basename;
use strict;

die "usage: $0 bgfFile cerius_ff [save name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($bgfFile, $ceriusFF, $saveName) = @ARGV;

FileTester($bgfFile);
FileTester($ceriusFF);

if (! $saveName) {
    $saveName = basename($bgfFile);
    $saveName =~ s/\.\w+$//;
    $saveName .= "_mod.bgf";
}

my ($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
my ($PARMS)= parseCerius2FF($ceriusFF);
my ($BOX) = GetBox($ATOMS, undef, undef);

my ($bInfo) = "PERIOD 111\nAXES   ZYX\nSGNAME P 1                  1    1\n" .
    sprintf("CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n", $BOX->{"X"}{"hi"} - $BOX->{"X"}{"lo"},
	    $BOX->{"Y"}{"hi"} - $BOX->{"Y"}{"lo"}, $BOX->{"Z"}{"hi"} - $BOX->{"Z"}{"lo"}, 
	    $BOX->{"X"}{"angle"}, $BOX->{"Y"}{"angle"}, $BOX->{"Z"}{"angle"}) . 
    "CELLS    -1    1   -1    1   -1    1";
push @{ $HEADERS }, $bInfo;
addHeader(\%{ $ATOMS }, $HEADERS);
createBGF($ATOMS, $BONDS, $saveName);
		     
