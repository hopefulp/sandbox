#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsTrjType GetLammpsByteOffset);
use Packages::General qw(FileTester TrjSelections GetSelections);
use Packages::AMBER qw(CreateAmberTrj);
use Packages::ManipAtoms qw(GetAtmList GetMols);
use File::Basename;
use Getopt::Std qw(getopt);

my ($velFile, $saveName, $saveType, $SELECT, $pStr, $bgfFile);
my ($ATOMS, $BONDS, $SOLVENT, $velFactor, $isShake, $ATOMSELECT);

sub init;
sub saveVel;
sub getShakeDOF;

$|++;
&init;
$velFactor = 1;
if (defined($bgfFile) and -e $bgfFile and -r $bgfFile and -T $bgfFile) {
    print "Parsing bgf file $bgfFile...";
    ($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
    $velFactor = getShakeDOF($ATOMS,$BONDS) if ($isShake);
    print "Done\n";
    if ($ATOMSELECT) {
	print "Parsing atom/residue selection...";
        $ATOMSELECT = GetSelections($ATOMSELECT, 0);
	$ATOMSELECT = GetAtmList($ATOMSELECT, $ATOMS);
	print "Done\n";
    } else {
	for (sort {$a<=>$b} keys %{ $ATOMS }) {
	    $ATOMSELECT->{$_} = $_;
	}
    }
}

$velFactor *= 1000/20.455 if ($saveType ne "lammps");

open OUTDATA, "> $saveName" || die "ERROR: Cannot create $saveName: $!\n";
$pStr = "Parsing LAMMPS Velocity file $velFile...";
if ($saveType eq "lammps") {
    print OUTDATA "Velocities\n\n";
} else {
    print OUTDATA "AMBER velocity file created by $ENV{USER} at " . 
	scalar(localtime) . "\n";
}
&GetLammpsByteOffset($SELECT, $velFile, scalar keys %{ $ATOMS });
&GetLammpsTrjType($SELECT, $velFile, undef);
&ParseLAMMPSTrj($saveType, $velFile, $SELECT, "vel", \&saveVel, $pStr, \*OUTDATA);
close OUTDATA or die "ERROR: Cannot close $saveName: $!\n";

sub getShakeDOF {
    my ($atoms, $bonds) = @_;
    my ($i, $dof, $factor, $atmCount);

    &GetMols($atoms, $bonds);
    print "adjusting velocities for dof removed by shake...";
    $atmCount = scalar keys %{ $atoms };
    $dof = $atmCount * 3;

    for $i (keys %{ $atoms }) {
	if ($atoms->{$i}{FFTYPE} =~ /^H/) { #if this is a hydrogen
	    $dof--;
	    if ($atoms->{$i}{MOLECULE}{MOLSIZE}== 3) {
		$dof -= 0.5;
	    }
	}
    }

    $factor = sqrt($atmCount * 3 / $dof);

    return $factor;
}

sub saveVel {
    my ($DATA, $type, $fileHandle) = @_;
    my ($i, $j, $BOX, $count);

    #%{ $ATOMSELECT } = (keys %{ $DATA->{ATOMS} }) if (! defined($ATOMSELECT));
    $count = 1;
    if ($type eq "lammps") {
	for $i (sort {$a<=>$b} keys %{ $ATOMSELECT }) {
	    next if (! exists($DATA->{ATOMS}{$i}));
	    printf $fileHandle "%8d", $i;
	    for $j ("XVEL", "YVEL", "ZVEL") {
		printf $fileHandle "%12.8f", $DATA->{ATOMS}{$i}{$j} * $velFactor;
	    }
	    print $fileHandle "\n";
	    $count++;
	}
    } else {
	for $i (1 .. $DATA->{"NUMBER OF ATOMS"}[0]) {
	    for $j ("X", "Y", "Z") {
		$DATA->{ATOMS}{$i}{$j . "COORD"} = $DATA->{ATOMS}{$i}{$j . "VEL"} * $velFactor;
	    }
	}
	CreateAmberTrj($DATA->{ATOMS},$BOX, $fileHandle);
    }
}

sub init {
    my (%OPTS, $selections, $atomSelect);
    
    getopt('vtsfbia',\%OPTS);
    ($velFile, $selections, $saveType, $saveName, $bgfFile, $isShake, $atomSelect) = 
	($OPTS{v}, $OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{b}, $OPTS{i}, $OPTS{a});
    die "usage:$0 -v lammps velocity file -f trajectory frames (* for all)\n -b [bgf file] -a (atom selection)" . 
	"\t-t ouput type [=lammps|amber] (optional) " . "-s [saveName] (optional) -i (is shake = no)\n"
	if (! defined($velFile) || ! defined($selections));

    print "Initializing...";
    FileTester($velFile);
    if (! defined($saveName)) {
	$saveName = basename($velFile);
	$saveName =~ s/\.\w+$/_lammps\.vel/;
    }
    $saveType = "lammps" if (! defined($saveType) || $saveType !~ /^[lammps|amber]/i);
    $saveType = lc($saveType);
    $SELECT = TrjSelections($selections);
    if (defined($atomSelect)) {
	if ($atomSelect =~ /\s+/) {
	    @{ $ATOMSELECT } = split /\s+/, $atomSelect;
	} else {
	    $ATOMSELECT->[0] = $atomSelect;
	}
    }

    $isShake = 0 if (! defined($isShake) or $isShake !~ /^(1|yes)$/i);
    print "Done\n";
}
