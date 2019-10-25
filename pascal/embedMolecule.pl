#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::General qw(FileTester CoM CombineMols LoadFFs);
use Packages::FileFormats qw(GetBGFFileInfo createBGF createHeaders addHeader addBoxToHeader GetBGFAtoms);
use Packages::BOX qw(GetBox CreateGrid GetNeighbours);
use Packages::ManipAtoms qw(GetMols MapOrigin SplitAtomsByMol);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub checkAtomTypes;
sub embeddMols;
sub removeOverlaps;
sub addRadii;
sub removeOverlapsN;

my ($DATA, $FFILES, $ATOMS, $BONDS, $HEADERS, $BOX, $PARMS);
my ($i, $soluBGF, $centerType, $solvBGF, $saveName, $solvMols);
my ($GRID, $BBOX, $atm_counter, $TMP, $soluAtms);

$|++;
$FFILES = &init;
$PARMS = LoadFFs($FFILES) if (defined($FFILES));
print "Parsing solvent/membrane bgf $solvBGF...";
($DATA->{SOLVENT}{ATOMS}, $DATA->{SOLVENT}{BONDS}, $DATA->{SOLVENT}{HEADERS}) = GetBGFFileInfo($solvBGF, 1);
print "Done\nParsing solute bgf $soluBGF...";
($DATA->{SOLUTE}{ATOMS}, $DATA->{SOLUTE}{BONDS}, $DATA->{SOLUTE}{HEADERS}) = GetBGFFileInfo($soluBGF, 1);
print "Done\n";
&checkAtomTypes($DATA, $PARMS) if (defined($FFILES));
for $i ("SOLUTE", "SOLVENT") {
    print "Computing ${centerType} for $i...\r";
    $DATA->{$i}{BOX} = GetBox($DATA->{$i}{ATOMS});
}
print "Computing ${centerType}s...Completed\n";
print "Embedding systems...";
&embeddMols($DATA) if ($centerType ne "NONE");
($ATOMS, $BONDS) = CombineMols($DATA->{SOLUTE}{ATOMS}, $DATA->{SOLVENT}{ATOMS},
			       $DATA->{SOLUTE}{BONDS}, $DATA->{SOLVENT}{BONDS});
&GetMols($ATOMS, $BONDS);
$BOX = GetBox($ATOMS, undef, $DATA->{SOLVENT}{HEADERS});
#&MapOrigin($BOX, CoM($ATOMS));
&addRadii($ATOMS, $PARMS);
print "Done\nCreating grid...\n";
($GRID, $BBOX, $atm_counter) = CreateGrid($ATOMS, 0, $BOX, 1, 1);
print "Removing overlaps...";
$soluAtms = scalar keys %{ $DATA->{SOLUTE}{ATOMS} };
$solvMols = SplitAtomsByMol($ATOMS, $ATOMS);
$soluAtms = 0;
$TMP = &removeOverlaps($ATOMS, $GRID, $soluAtms, $solvMols);
($ATOMS, $BONDS) = GetBGFAtoms($TMP, $ATOMS, $BONDS);
print "Done\nCreating BGF file $saveName...";
$HEADERS = createHeaders(undef, $saveName);
&addBoxToHeader($HEADERS, $BOX);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub removeOverlaps {
    my ($bgfAtoms, $grid, $offset, $molList) = @_;
    my ($i, $x, $y, $z, $CELL, $solvent, $CLIST, $molID); 
    my ($dist, $index, $currCell, $atoms, $j, $solAtom);

    %{ $atoms } = %{ $bgfAtoms };
    for $i (keys %{ $atoms }) {
 	next if (! exists($atoms->{$i}) or $atoms->{$i}{IS_SOLVENT});
	$x = $atoms->{$i}{CELL}{"XINDEX"};
	$y = $atoms->{$i}{CELL}{"YINDEX"};
	$z = $atoms->{$i}{CELL}{"ZINDEX"};
	next if (! $x or ! $y or ! $z);
	$currCell = $grid->{$x}{$y}{$z};
	next if (exists($currCell->{VISITED}));
	$CLIST = GetNeighbours($GRID, $currCell);
	$currCell->{VISITED} = 1;
	for $solAtom (@{ $currCell->{ATOMS} }) {
	    for $CELL (@{ $CLIST }) {
		for $solvent (@{ $CELL->{WATERS} }) {
		    next if (! $solvent->{IS_SOLVENT});
		    $dist = 0;
		    for $index ("XCOORD", "YCOORD", "ZCOORD") {
			$dist += (($solAtom->{$index} - $solvent->{$index}) ** 2);
		    }
		    $dist = sqrt($dist);
		    if ($dist < 0.5) { #remove entire solvent molecule
			$molID = $solvent->{MOLECULEID};
			for $j (keys %{ $molList->{$molID} }) {
			    delete $atoms->{($j + $offset)};
			}
		    }
		}
	    }
	}
    }
    
    return $atoms;
}

sub init {
    my (%OPTS, @tmp, @FFILES, $i, $FF);
    
    getopt('msfcw',\%OPTS);
    ($solvBGF, $soluBGF, $FF, $centerType, $saveName) = ($OPTS{m}, $OPTS{s}, $OPTS{f}, $OPTS{c}, $OPTS{w});
    
    die "usage: $0 -m solvent/membrane bgf -s solute bgf -f [\"forcefield1 forcefield2...\"]\n" . 
	"\t-c [com|cog|none] centering type -w [saveName.bgf]\n"
	if (! defined($solvBGF) || ! defined($soluBGF));
    
    print "Initializing...";
    FileTester($solvBGF);
    FileTester($soluBGF);

    $centerType = "COM" if (! defined($centerType) or $centerType !~ /^com|cog|cop|none/);
    $centerType = uc($centerType);
    $centerType = "COP" if ($centerType !~/COM|NONE/);

    if (! defined($saveName)) {
	$saveName = basename($soluBGF);
	$saveName =~ s/\.\w+$/_embed\.bgf/;
    }

    if ($FF =~ /\s+/) {
        @tmp = split /\s+/, $FF;
    } else {
        @tmp = ($FF);
    }
  
    for $i (@tmp) {
        if (-e $i &&  -r $i && -T $i) {
            push @FFILES, $i;
	}
    }

    die "ERROR: No valid CERIUS2 forcefield found!\n" if (! @FFILES);
    print "Done\n";
    return \@FFILES;
}

sub addRadii {
    my ($atoms, $parms) = @_;
    my ($i, $ffType);

    for $i (keys %{ $atoms }) {
	$ffType = $atoms->{$i}{FFTYPE};
	if (! defined($PARMS)) {
	    $atoms->{$i}{RADII} = 1;
	} else {
            $atoms->{$i}{RADII} = $parms->{VDW}{$ffType}{$ffType}{1}{VALS}[1];
	}
    }
}

sub checkAtomTypes {
    my ($mols, $parms) = @_;
    my ($i, $j, $ffType);

    for $i ("SOLUTE", "SOLVENT") {
	for $j (keys %{ $mols->{$i}{ATOMS} }) {
	    $ffType = $mols->{$i}{ATOMS}{$j}{FFTYPE};
	    die "ERROR: Force field type $ffType not found in forcefield(s). Aborting\n" 
		if (! exists($parms->{ATOMTYPES}{$ffType}));
	    if ($centerType eq "COM") {
		$mols->{$i}{ATOMS}{$j}{RADII} = $parms->{VDW}{$ffType}{$ffType}{1}{VALS}[1];
		$mols->{$i}{ATOMS}{$j}{MASS} = $parms->{ATOMTYPES}{$ffType}{MASS};
	    }
	    if ($i eq "SOLVENT") {
		$mols->{$i}{ATOMS}{$j}{IS_SOLVENT} = 1;
	    } else {
		$mols->{$i}{ATOMS}{$j}{IS_SOLUTE}  = 1;
	    }
	}
    }
}

sub embeddMols_old {
    my ($mols) = $_[0];
    my ($i, $j, @dims);

    @dims = keys %{ $mols->{SOLVENT}{BOX} };
    for $i (@dims) {
	$mols->{SOLUTE}{OFFSET}{$i} = (($mols->{SOLUTE}{BOX}{$i}{hi} + $mols->{SOLUTE}{BOX}{$i}{lo})/2 - 
				      ($mols->{SOLVENT}{BOX}{$i}{hi} + $mols->{SOLVENT}{BOX}{$i}{lo})/2);
    }

    for $i (keys %{ $mols->{SOLUTE}{ATOMS} }) {
	for $j (@dims) {
	    $mols->{SOLUTE}{ATOMS}{$i}{"${j}COORD"} -= $mols->{SOLUTE}{OFFSET}{$j};
	}
    }
}

sub embeddMols {
    my ($mols) = $_[0];
    my ($i, $j, @dims, $coms);

    $coms->{SOLUTE} = CoM($mols->{SOLUTE}{ATOMS});
    $coms->{SOLVENT} = CoM($mols->{SOLVENT}{ATOMS});

    @dims = keys %{ $mols->{SOLVENT}{BOX} };
    for $i (@dims) {
	$mols->{SOLUTE}{OFFSET}{$i} = ($coms->{SOLUTE}{"${i}COORD"} - $coms->{SOLVENT}{"${i}COORD"});
    }

    for $i (keys %{ $mols->{SOLUTE}{ATOMS} }) {
	for $j (@dims) {
	    $mols->{SOLUTE}{ATOMS}{$i}{"${j}COORD"} -= $mols->{SOLUTE}{OFFSET}{$j};
	}
    }
}
