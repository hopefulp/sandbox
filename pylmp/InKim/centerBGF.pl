#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF AddMass);
use Packages::General qw(FileTester CoM LoadFFs);
use Packages::BOX qw(GetBox MoveAtomsToOrigin);
use Getopt::Std qw(getopt);
use Packages::ManipAtoms qw(CenterSystem);

sub init;
sub showUsage;
sub getCenter;

my ($bgfFile, $cerius2FF, $saveName, $centerType);
my (%Offset, $dim, $atomC, $radii, $CENTER);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
my ($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nParsing CERIUS2 forcefield $cerius2FF...";
my ($PARMS) = LoadFFs($cerius2FF);
AddMass($ATOMS, $PARMS);
print "Done\nGetting box information...";
my ($BOX) = GetBox($ATOMS, $PARMS, $HEADERS);
$CENTER = getCenter($ATOMS, $BOX, $centerType);
print "Done\nCentering atoms to origin...";
if ($centerType eq "box_origin") {
    &MoveAtomsToOrigin($ATOMS);
} else {
    CenterSystem($ATOMS, $CENTER, ());
}
addHeader($ATOMS, $HEADERS);
print "Done\nCreating BGF file $saveName...";
createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub getCenter {
    my ($atoms, $box, $centerOpt) = @_;
    my ($atomCenter) = CoM($atoms);
    my (%CENTER, $i);
    my (@dim) = ("X", "Y", "Z");

    if ($centerOpt eq "box_origin") {
	$box = GetBox($ATOMS, $PARMS, undef);
	for $i (@dim) {
	    $CENTER{"${i}COORD"} = $box->{$i}{lo};
	}
    } elsif ($centerOpt eq "com_origin") {
	for $i (@dim) {
	    $CENTER{"${i}COORD"} = $atomCenter->{"${i}COORD"} - $box->{$i}{lo};
	}
    } elsif ($centerOpt eq "box_center") {
	for $i (@dim) {
	     $CENTER{"${i}COORD"} = -1*$box->{$i}{len}/2 - $box->{$i}{lo};
	}
    } elsif ($centerOpt eq "com_center") {
	for $i (@dim) {
	     $CENTER{"${i}COORD"} = -($box->{$i}{len}/2 - $box->{$i}{lo} - $atomCenter->{"${i}COORD"});
	}
    } else {
	for $i (@dim) {
	    $CENTER{"${i}COORD"} = $box->{$i}{len}/2 - $atomCenter->{"${i}COORD"} + $box->{$i}{lo};
	    #$CENTER{"${i}COORD"} = $atomCenter->{"${i}COORD"};
	}
    }

    return \%CENTER;
}

sub init {
    my (%OPTS, $usage, $centerOpt, $ffData);
    getopt('bfsc',\%OPTS);
    ($bgfFile, $ffData, $saveName, $centerOpt) = ($OPTS{b},$OPTS{f},$OPTS{s},$OPTS{c});
    $usage = &showUsage;
    for ($bgfFile, $ffData) {
	die "$usage\n" if (! defined($_));
    }

    print "Initializing...";
    FileTester($bgfFile);
    while ($ffData =~ /(\S+)/g) {
	if (-e $1 and -r $1 and -T $1) {
	    push @{ $cerius2FF }, $1;
	}
    }
    die "ERROR: No valid force found while search \"$ffData\"\n" if (! $cerius2FF);

    $centerOpt = "box_origin" if (! defined($centerOpt));
    
    if ($centerOpt =~ /(box_origin|com_center|com_origin|box_center)/i) {
	$centerType = lc($1);
    }
    if (! $saveName) {
	$saveName = $bgfFile;
	$saveName =~ s/\.\w+$/_${centerType}\.bgf/;
    }

    die "ERROR: Invalid centering option \"${centerOpt}\"\n$usage\n" if (! defined($centerType));
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -f cerius2 forcefield -s [save name] -c [centering type]\n" .
	"options:\n\t-b bgf file: location of the bgf file\n" .
	"\t-f cerius2 forcefield: location of CERIUS2 formatted forcefield.\n\t\tBGF file must be typed with this forcefield\n" .
	"\t-s [save name]: (optional) the name to save the centered bgf. \n\t\tIf not specified will be {bgf file}_centered.bgf\n" .
	"\t-c [centering type]: (optional) the type of centering to perform.\n" .
	"\t\tChoices:\n\t\tcom_origin - center the center of mass of the molecule to the origin\n" .
	"\t\tcom_center - place the center of mass of the molecule in the center of the box\n" .
	"\t\tbox_origin - (default) place the lowest extrema of the molecule at the origin\n";
    return $usage;
}
