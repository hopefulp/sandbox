#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::General qw(CoM FileTester GetSelections LoadFFs);
use Packages::FileFormats qw(GetBGFFileInfo AddMass);
use Packages::CERIUS2 qw(parseCerius2FF);
use Getopt::Std qw(getopt);
use Packages::ManipAtoms qw(GetAtmList GetAtmData);

sub init;

my ($bgfFile, $cerius2FF, $SELECT);
my ($ATOMS, $BONDS, $PARMS, $CENTER);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
print "Done\n";
if (defined($cerius2FF)) {
    $PARMS = LoadFFs($cerius2FF);
    AddMass($ATOMS, $PARMS);
}
print "Center of Mass: ";
$CENTER = CoM(GetAtmData($ATOMS, GetAtmList($SELECT, $ATOMS)));
#$CENTER = CoM($ATOMS);
for ("XCOORD", "YCOORD", "ZCOORD") {
    printf "%15.5f", $CENTER->{$_};
}
print "\n";

sub init {
    my (%OPTS, $select, $selection, $fflist);
    
    getopt('bfa',\%OPTS);
    ($bgfFile, $fflist, $select) = ($OPTS{b},$OPTS{f}, $OPTS{a});
    
    die "usage: $0 -b bgf file -f cerius2 force field -a (atom selection. default all)\n" 
	if (! defined($bgfFile));

    print "Initializing...";
    FileTester($bgfFile);
    if (defined($fflist)) {
        while ($fflist =~ /(\S+)/g) {
	    if (-e $1 and -r $1 and -T $1) {
		push @{ $cerius2FF }, $1;
	    }
	}
	die "ERROR: No valid forcefield files found while search \"$fflist\"\n" if (! $cerius2FF);
    }
    $select = "*" if (! defined($select));
    if ($select =~ /\s+/) {
	@{ $selection } = split /\s+/, $select;
    } else {
	$selection->[0] = $select;
    }

    $SELECT = GetSelections($selection, 0);

}
