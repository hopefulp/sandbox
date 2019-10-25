#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use Packages::General qw(FileTester TrjSelections Trim CenterOnMol CoM GetSoluteAtoms GetBondLength);
use Packages::FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use Packages::LAMMPS qw(ParseLAMMPSTrj CreateLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use Packages::AMBER qw(CreateAmberTrj);
use Packages::ManipAtoms qw(UnwrapAtoms ScaleAtoms ImageAtoms SelectAtoms BuildAtomSelectionString  
				SplitAtomsByMol GetAtmData GetMols);
use Packages::BOX qw(MoveAtomsToOrigin);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub saveCoords;
sub init;
sub createLAMMPSHeader;
sub updateBGF;
sub doResFix;
sub showUsage;
sub createCOMTrj;
sub createCOMBGF;
sub createCOMAtom;
sub writeCOMBGF;
sub getSoluteSolvent;
sub intAMBERvels;

my ($lammpsTrj, $bgfFile, $selection, $reImageCoord, $saveName, $outType, $comTrj, $slabDim);
my ($printAtoms, $SELECT, $LAMMPSOPTS, $pStr, $SOLUTEATMS, $count, $reScaleCoord, $MOLS);
my ($BGF, $BONDS, $soluteSelect, $isTrj, $molList, $comBGF, $imgPoint, $savePt, $com_center);
my ($AMBERVELHANDLE, $writeAMBERvels);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($BGF, $BONDS) = GetBGFFileInfo($bgfFile, 0);
$MOLS = GetMols($BGF, $BONDS);
($SOLUTEATMS, $molList) = getSoluteSolvent($BGF, $BONDS, $soluteSelect, $MOLS);
#die "ERROR: No valid solute atoms found!\n" if (! keys %{ $SOLUTEATMS });
#die "ERROR: No valid solvent atoms found!\n" if (! keys %{ $molList });
$comBGF = createCOMBGF($BGF) if ($comTrj);
print "Done\n";
&GetLammpsByteOffset($SELECT, $lammpsTrj, scalar keys %{ $BGF });
&GetLammpsTrjType($SELECT, $lammpsTrj, "atom", \%{ $LAMMPSOPTS });
$AMBERVELHANDLE = &initAMBERvels($saveName) if ($isTrj and $outType eq "AMBER" and $LAMMPSOPTS->{hasvel});
$pStr = "Parsing LAMMPS trajectory $lammpsTrj...";
$count = 0;
&ParseLAMMPSTrj($BGF, $lammpsTrj, $SELECT, "atom", \&saveCoords, $pStr, \*OUTDATA);
close OUTDATA or die "ERROR: Cannot close $saveName: $!\n" if ($isTrj);
close $AMBERVELHANDLE or die "ERROR: Cannot close AMBER velocity file: $!\n" if ($isTrj and $writeAMBERvels);
&writeCOMBGF($comBGF, $saveName) if ($comTrj);

sub initAMBERvels {
    my ($coordname) = @_;
    my ($velname) = $coordname;
    $velname =~ s/\.(\w)+$//;
    $velname .= ".mdvel";
    $velname = ".vel.mdcrd" if($1 eq "mdvel");

    open AMBERVEL, "> $velname" or die "ERROR: Cannot create AMBER velocity file $velname: $!\n";
    print AMBERVEL "TITLE: $outType Velocity file create by lammpsTrj2amber on " . scalar(localtime(time())) . "\n";

    $writeAMBERvels = 1;
    return \*AMBERVEL;
}

sub saveCoords {
    my ($DATA, $bgfInfo, $fileHandle) = @_;
    my ($BOX, $mID, $HEADERS, @tmp, $CENTER, $j, $tstep); 
    my ($MOLECULE, $i, $fName, $index, $totAtms);

    $tstep = $DATA->{"TIMESTEP"}[0];
    $BOX = ConvertLammpsBox($DATA->{"BOX BOUNDS"});
    $totAtms = $DATA->{"NUMBER OF ATOMS"}[0];
    if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
	UnwrapAtoms($DATA->{ATOMS}, $BOX, $LAMMPSOPTS->{scaled});
    }
    if ($reImageCoord) {
        @tmp = ("XCOORD", "YCOORD", "ZCOORD");
        for $j (@tmp) {
            $CENTER->{$j} = $BOX->{$j}{CENTER} = $BOX->{$j}{len}/2;
	}
        if (defined($SOLUTEATMS)) {
	    $MOLECULE = GetAtmData($DATA->{ATOMS}, $SOLUTEATMS);
            $CENTER = CoM($MOLECULE);
	} elsif (defined($imgPoint)) {
	    $CENTER = $imgPoint;
	} elsif (defined($com_center)) {
	    $CENTER = CoM($DATA->{ATOMS});
	}
	for $j (@tmp) {
	    $BOX->{$j}{hi} = $BOX->{$j}{len};
	    $BOX->{$j}{lo} = 0;
	    $CENTER->{$j} = $BOX->{$j}{CENTER} if(defined($slabDim) and $j eq $slabDim);
	    for $i (keys %{ $DATA->{ATOMS} }) {
		$DATA->{ATOMS}{$i}{$j} += ($BOX->{$j}{CENTER} - $CENTER->{$j});
	    }
	}
        for $index (keys %{ $MOLS }) {
	    $MOLECULE = GetAtmData($DATA->{ATOMS}, $MOLS->{$index}{MEMBERS});
	    $CENTER = CoM($MOLECULE);
	    ImageAtoms($MOLECULE, $CENTER, $BOX);
        }
	if ($savePt) {
	    %{ $DATA->{ATOMS}{$totAtms+1} } = %{ $imgPoint };
	}
    }
    $DATA->{ATOMS} = createCOMTrj($DATA->{ATOMS}, $comBGF) if ($comTrj);
    if ($outType eq "AMBER") {
	&CreateAmberTrj($DATA->{ATOMS},$BOX, $fileHandle);
	if($writeAMBERvels) {
	   for $i (keys %{ $DATA->{ATOMS} }) {
	     for $j ("X", "Y", "Z") {
		$DATA->{ATOMS}{$i}{"${j}COORD"} = $DATA->{ATOMS}{$i}{"${j}VEL"} * 1000/20.455; #copy vel info into coord field
	     }
	   }
	   &CreateAmberTrj($DATA->{ATOMS},undef, $AMBERVELHANDLE);
	}
    } elsif ($outType eq "LAMMPS") {
	ScaleAtoms($DATA->{ATOMS}, $BOX) if ($reScaleCoord);
	#$HEADERS = createLAMMPSHeader($DATA, $BOX);
        if ($reScaleCoord) {
	    push @{ $HEADERS->{cType} }, "xs ys zs";
	} else {
	    push @{ $HEADERS->{cType} }, "xu yu zu";
	}
	CreateLAMMPSTrj($DATA->{ATOMS}, $DATA, $fileHandle);
    } else { #bgf
	$count++;
	$fName = $saveName;
	@tmp = keys %{ $SELECT };
	if (defined($SELECT) && $#tmp > 0) {
	    $fName =~ s/\.\w+$//;
	    $fName .= ".${tstep}.bgf";
	}
	updateBGF($bgfInfo, $DATA->{ATOMS}, $BOX, $DATA->{TIMESTEP}[0]);
	createBGF($bgfInfo, $BONDS, $fName);
    }
}

sub createCOMAtom {
    my ($index) = $_[0];
    my (%ATOM);

    $ATOM{"INDEX"} = $index;
    $ATOM{"ATMNAME"} = "Na";
    $ATOM{"RESNAME"} = "Na";
    $ATOM{"RESNUM"} = $index;
    $ATOM{"FFTYPE"} = "Na";
    $ATOM{"NUMBONDS"} = 0;
    $ATOM{"LONEPAIRS"} = 0;
    $ATOM{"CHARGE"} = 0;
    $ATOM{"RADII"} = 0;
    $ATOM{"LABEL"} = "HETATM";
    return \%ATOM;
}

sub createCOMBGF {
    my ($atoms) = $_[0];
    my ($mols, $i, $comAtoms, $comBGF, $j, $atomCoord);

    $mols = SplitAtomsByMol($BGF, $BGF);

    for $i (keys %{ $mols }) {
        $comBGF->{$i} = createCOMAtom($i);
	$comAtoms = GetAtmData($atoms, $mols->{$i});
	for $j (keys %{ $mols->{$i} }) {
	    $comBGF->{$i}{MOLECULE}{MEMBERS}{$j} = $j;
	}
	$atomCoord = CoM($comAtoms);
	for $j (keys %{ $atomCoord }) {
	    $comBGF->{$i}{$j} = $atomCoord->{$j};
	}
    }

    return $comBGF;
}

sub createCOMTrj {
    my ($atomCoords, $comAtoms) = @_;
    my ($i, $atomCoord, $j, $atoms);

    for $i (keys %{ $comAtoms }) {
	$atoms = GetAtmData($atomCoords, $comAtoms->{$i}{MOLECULE}{MEMBERS});
	$atomCoord = CoM($atoms);
	$comAtoms->{$i}{TYPE} = 1;
	for $j (keys %{ $atomCoord }) {
	    $comAtoms->{$i}{$j} = $atomCoord->{$j};
	}
    }

    return $comAtoms;
}

sub writeCOMBGF {
    my ($comAtoms, $trjSaveName) = @_;
    my ($bgfSaveName, $i, $comBonds, $comHeaders);

    $bgfSaveName = $trjSaveName;
    $bgfSaveName =~ s/\.\w+$//;
    $bgfSaveName .= "_com.bgf";
    print "Creating COM BGF file $bgfSaveName...";
    for $i (keys %{ $comAtoms }) {
	$comBonds->{$i} = ();
    }
    $comHeaders = createHeaders(undef, $bgfSaveName);
    &addHeader($comAtoms, $comHeaders);
    &createBGF($comAtoms, $comBonds, $bgfSaveName);
    print "Done\n";
}

sub doResFix {
    my ($ATOMS, $BOX, $currAtm) = @_;
    my ($i, @dims, $j, $bondPos, $atomPos);

    $ATOMS->{$currAtm}{VISITED} = 1;
    @dims = ("X", "Y", "Z");
    for $j (@dims) {
	$atomPos->{$j} = $ATOMS->{$currAtm}{$j . "COORD"};
    }
    for $i ( @{ $BONDS->{$currAtm} }) {
	next if (exists($ATOMS->{$i}{VISITED}));
	for $j (@dims) {
	    $bondPos = $ATOMS->{$i}{$j . "COORD"};
	    if(abs($bondPos - $atomPos->{$j}) > 5) { #image problem
		if ($atomPos->{$j} > $bondPos) { # move in + direction
		    while (($atomPos->{$j} - $bondPos) > 5) {
			$bondPos += ($BOX->{$j}{hi} - $BOX->{$j}{lo});
		    }
		} else { # move in - direction
                    while (($bondPos - $atomPos->{$j}) > 5) {
                        $bondPos -= ($BOX->{$j}{hi} - $BOX->{$j}{lo});
                    }
		}
		$ATOMS->{$i}{$j . "COORD"} = $bondPos;
		$ATOMS->{$i}{$j . "INDEX"} = $ATOMS->{$currAtm}{$j . "INDEX"};
	    }
	}
	doResFix($ATOMS, $BOX, $i);
    }
}	     
	    
sub updateBGF {
    my ($REF, $MOD, $BOX, $TSTEP) = @_;
    my (@tmp, $i, $HEADERS, $j);

    $BOX->{1}{DATA} = 90;
    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (2 .. 4) {
	$BOX->{$i}{DATA} = ($BOX->{$tmp[($i - 2)]}{hi} - $BOX->{$tmp[($i - 2)]}{lo});
    }

    $TSTEP = "ts_${TSTEP}";
    $HEADERS = createHeaders($BOX, $TSTEP);

    for $i (keys %{ $REF }) {
	for $j (keys %{ $MOD->{$i} }) {
	    $REF->{$i}{$j} = $MOD->{$i}{$j};
	}
    }
    addHeader($REF, $HEADERS);
}

sub getSoluteSolvent {
    my ($atoms, $bonds, $selection, $mols) = @_;
    my ($solute, $solventMols, $i, $solventList);
    
    $solute = SelectAtoms($selection, $atoms) if (defined($selection));
    for $i (keys %{ $atoms }) {
	next if (defined($solute) && exists($solute->{$i}));
	$solventMols->{${ $atoms->{$i}{MOLECULEID} }} = $atoms->{$i}{MOLECULE};
    }
    return ($solute, $solventMols);
}


sub init {
    my ($i, %OPTS, $usage, $atomSelection, @tmp, $coordOpt, $pointStr);

    getopt('lstboamcpqdn',\%OPTS);
    $usage = &showUsage;

    ($lammpsTrj,$bgfFile,$selection,$outType,$saveName,$coordOpt,$atomSelection,$comTrj,$pointStr,$savePt,$slabDim, $com_center) = 
    ($OPTS{l}, $OPTS{b}, $OPTS{t}, $OPTS{o}, $OPTS{s}, $OPTS{m}, $OPTS{a}, $OPTS{c}, $OPTS{p}, $OPTS{q}, $OPTS{d}, $OPTS{n});

    for ($lammpsTrj, $bgfFile, $selection) {
	die "$usage\n" if (! defined($_));
    }
    print "Initializing...";

    FileTester($lammpsTrj);
    FileTester($bgfFile);
    
    $outType = "lammps" if (! defined($outType));
    if ($outType !~ /(lammps|amber|bgf)/i) {
	print "WARNING: Expected lammps, amber or bgf for out_type. Got $outType...\n" . 
	    "Setting to default value: lammps\n";
	$outType = "lammps";
    } else {
	$outType = $1;
    }
    print "$outType traj ";

    $coordOpt = 0 if (! defined($coordOpt) or $coordOpt !~ /^(1|2|3)$/);
    if ($coordOpt == 0) {
	$reScaleCoord = 0;
	$reImageCoord = 0;
    } elsif($coordOpt==1) {
        $reScaleCoord = 0;
        $reImageCoord = 1;
    } elsif($coordOpt==2) {
        $reScaleCoord = 1;
        $reImageCoord = 1;
    } else {
        $reScaleCoord = 1;
        $reImageCoord = 0;
    }
    print "with coordinate imaging.." if ($reImageCoord);
    print "without coordinate imaging..." if (! $reImageCoord);

    if (! defined($saveName)) {
	$saveName = basename($lammpsTrj);
	$saveName =~ s/\.\w+$/_mod\.lammpstrj/;
    }

    $SELECT = TrjSelections($selection);
    $outType = uc($outType);
    
    $isTrj = 0;
    $isTrj = 1 if ($outType =~ /LAMMPS|AMBER/i);
    if ($isTrj) {
	open OUTDATA, "> $saveName" or die "Cannot write to trajectory file $saveName: $!\n";
	$pStr = "Creating $outType trajectory $saveName...";
	print OUTDATA "TITLE: $outType Trajectory file create by lammpsTrj2amber on " . scalar(localtime(time())) . "\n";
    } else {
	$pStr = "Creating $outType files for each frame...";
    }

    if (defined($atomSelection)) {
	$soluteSelect = BuildAtomSelectionString($atomSelection);
    }
    if (defined($pointStr) and $pointStr =~ /(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)/) {
	($imgPoint->{XCOORD}, $imgPoint->{YCOORD}, $imgPoint->{ZCOORD}) = ($1, $2, $3);
	if (defined($savePt) && $savePt =~ /^(1|yes)/) {
	    $savePt = 1;
	} else {
	    $savePt = 0;
	}
    }
    undef($soluteSelect) if defined($imgPoint);
    $comTrj = 0 if (! defined($comTrj) or $comTrj !~ /^\s*(1|yes)\s*$/i);
    $comTrj = 1 if ($comTrj =~ /^\s*(1|yes)\s*$/i);
    if(defined($slabDim) && $slabDim =~ /^(x|y|z)$/) {
	$slabDim = uc $1 . "COORD";
    }

    $com_center = 0 if (! defined($com_center) or $com_center !~ /1|yes/i);
    $com_center = 1 if ($com_center =~ /1|yes/i);

    $writeAMBERvels = 0;
    print "Done\n";
}

sub showUsage {
    my ($usage) = <<DATA;
usage: $0 -b bgf file -l lammps trj -t trajectory selection -n [center com 2 box center = no] -m [reimage/rescale coordinates] -p [spacial point] -a [solute atoms] -o [output type] -s [save name] -c [create com trj] -q [save spacial point]
options:
	-b bgf file (Required): the location of the bgf file
	-l lammps trj (Required): the location of the lammps trajectory file
	-t trajectory selection: The number of frames to use. Can be a single integer, or several integers in quotes
		To specify a range specify it as :Ita-b:c, which will take frames a - b, every c. Specify multiple 
		ranges or a combination of ranges and single frames by enclosing them in quotes. \"*\" for all frames
	-m [reimage/rescale coordinates]: (Optional) Whether to reimage and rescale the coordinates back into the 
		unit cell. 0:No reimage or rescale (Default), 1:reimage no rescale, 2:rescale and rescale 
		3:rescale no reimage
	-p [spacial point]: (Optional) Image the coordinates wrt to point in space - or -
	-a [solute atoms]: (Optional) Image the coordinates wrt to the center of mass of certain atoms
		Any valid field (fftype, resname, xcoord etc) expression (e.g. "resname eq 'WAT'")
	-q [save spacial point]: (Optional) Save the coordinates of the spacial point. Defaul no
	-d [slab dimension]: (Optional) The non periodic dimension in a slab simulation
        -o [output type]: (Optional). Either lammps (default) or amber or bgf. If bgf, will print file for every frame

	-c [create com trj]: convert each molecule to an atom based on their center of mass. Default no
	-s [save name]: (Optional). The name the save the output trajectory
DATA

    return $usage;
}
