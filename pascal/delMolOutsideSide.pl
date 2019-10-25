#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use Packages::General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub removeMolsOutsideCell;
sub getCellBox;

my ($bgfFile, $saveFile);
my ($ATOMS, $BONDS, $HEADERS, $BOX);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$BOX = getCellBox($HEADERS);
print "Done\nRemoving any atoms outside unit cell...";
&removeMolsOutsideCell($ATOMS, $BONDS, $BOX);
print "Done\nCreating BGF file $saveFile...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub removeMolsOutsideCell {
    my ($atoms, $bonds, $boxInfo) = @_;
    my ($i, $dim, $j);

    for $i (keys %{ $atoms }) {
	next if (! exists($atoms->{$i}));
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    if ($atoms->{$i}{$dim} < 0 or $atoms->{$i}{$dim} > $boxInfo->{$dim}) {
		for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
		    delete $atoms->{$j};
		    delete $bonds->{$j};
		}
		delete $atoms->{$i} if (defined($atoms->{$i}));
		delete $bonds->{$i} if (defined($atoms->{$i}));
		last;
	    }
	}
    }
}

sub getCellBox {
    my ($headerInfo) = $_[0];
    my ($i, %CellBOX);

    for $i (@{ $headerInfo }) {
	if ($i =~ /^CRYSTX\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
	    $CellBOX{XCOORD} = $1;
	    $CellBOX{YCOORD} = $2;
	    $CellBOX{ZCOORD} = $3;
	}
    }

    die "ERROR: No Cell information found!\n" if (! %CellBOX);
    return \%CellBOX;
}

sub init {
    my (%OPTS);
    
    getopt('bs', \%OPTS);
    die "usage: $0 -b bgf file -s (save name)\n" if (! defined($OPTS{b}));
    
    print "Initializing...";
    ($bgfFile, $saveFile) = ($OPTS{b}, $OPTS{s});
    FileTester($bgfFile);
    
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_mod.bgf";
    }
    print "Done\n";
}
