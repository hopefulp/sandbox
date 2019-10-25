#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;

use Storable qw(dclone);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

use Packages::FileFormats qw(GetBGFFileInfo createBGF addHeader createPQR PrintCoord
			     DeleteAtoms InsertMol PrintCoord GetSystemCharge GetBGFAtoms);
use Packages::General qw(FileTester LoadFFs GetSelections IsDecimal CoM);
use Packages::CERIUS2 qw(ReadFFs);
use Packages::ManipAtoms qw(SplitAtomsByMol GetAtmList GetMols GroupAtomsByField 
			     CenterSystem MakeFieldSequential);
use Packages::BOX qw(GetBox);

use constant TOLERANCE => 0.0001;

sub init;
sub checkAtomTypes;
sub addRadii;
sub updateSolvent;
sub usage;
sub getIonParms;
sub getIonOpts;
sub getIonType;
sub determineIonPlacement;
sub updateIonFields;
sub placeIon;
sub numerically { ($a<=>$b); }
sub updateIonCoords;
sub calcElecPot;

my ($bgfFile, $saveName, $ffFiles, $solvSelect, $ionType, $ionConc);
my ($ATOMS, $BONDS, $HEADERS, $PARMS, $SOLVENT, $IONS, $BOX);
my ($soluCharge, $solvEng, $sMOLS, $numSolv);

$|++;
my ($start) = time();
$IONS = &init;
$PARMS = LoadFFs($ffFiles);
print "Getting Ion parameters...";
&getIonParms($IONS, $PARMS);
print "Done\nParsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
&GetMols($ATOMS,$BONDS);
print "Done\nGetting system dimensions...";
$BOX = GetBox($ATOMS, $PARMS, $HEADERS);
print "Done\nObtaining Solvent atoms...";
&checkAtomTypes($ATOMS, $PARMS);
($SOLVENT, undef) = GetAtmList($solvSelect, $ATOMS);
$sMOLS = GetMols($ATOMS, $BONDS, $SOLVENT);
$soluCharge = &updateSolvent($ATOMS, $SOLVENT);
$numSolv = scalar(keys %{ $sMOLS });
print "Found " . scalar(keys %{ $SOLVENT }) . " atoms in $numSolv molecules..Done\n";
&getIonOpts($IONS, $BOX, $ionConc, $soluCharge, $numSolv);
print "Calculating Electrostatic Potential...";
$solvEng = calcElecPot($ATOMS, $BONDS, $BOX, $HEADERS);
print "Done\nDetermining Ion placement based on potential...";
&determineIonPlacement($IONS, $solvEng);
print "Done\n";
&placeIons($ATOMS, $BONDS, \%{ $IONS->{1} }, $sMOLS);
if (exists($IONS->{2})) {
    &placeIons($ATOMS, $BONDS, \%{ $IONS->{2} }, $sMOLS);
}
print "Creating BGF file $saveName...";
&GroupAtomsByField($ATOMS,$BONDS,"RESNUM");
&MakeFieldSequential($ATOMS,"RESNUM");
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

$soluCharge = GetSystemCharge($ATOMS);
printf "Total System Charge: %.3f\n", $soluCharge;

my ($end) = time();

print "Elapsed: " . ($end - $start) . " seconds\n";

sub placeIons {
    my ($atoms, $bonds, $ions, $solvMols) = @_;
    my ($currMol, $solvAtoms, $currAtom, $i, $index, $rindex, $count); 
    my ($resId, $resName, $currIon, $resname, $soluCoM, @tmp);

    $index = 0;
    $rindex = 99999999;
    $count = $#{ $ions->{list} } + 2;
    for $i (keys %{ $atoms }) {
	$index = $i if ($i > $index);
	$rindex = $atoms->{$i}{RESNUM} if ($atoms->{$i}{IS_SOLVENT} and $atoms->{$i}{RESNUM} < $rindex);
	$atoms->{$i}{RESNUM} += $count if ($atoms->{$i}{IS_SOLVENT});
    }

    for $currMol (@{ $ions->{list} }) {
	$solvAtoms = ();
	@tmp = keys %{ $solvMols->{$currMol}{MEMBERS} };
	$currAtom = shift @tmp;
	next if (! exists($atoms->{$currAtom}) || ! exists($atoms->{$currAtom}{MOLECULE}{MEMBERS}));
	for $i (keys %{ $atoms->{$currAtom}{MOLECULE}{MEMBERS} }) {
	    $solvAtoms->{$i} = $atoms->{$i};
	}
        $currAtom = $atoms->{$currAtom};
        $resId = $currAtom->{"RESNUM"};
        $resname = $currAtom->{"RESNAME"};
        $rindex++;
        $soluCoM = CoM($solvAtoms);
	if(! exists($ions->{custom})) {
	    $currIon = updateIonFields($ions, $soluCoM, $index, $rindex);
	} else {
	    $currIon = updateIonCoords($ions, $soluCoM, $index, $rindex);
	}
        &DeleteAtoms($solvAtoms, $atoms, $bonds);
	&InsertMol($currIon->{atoms}, $currIon->{bonds}, $atoms, $bonds, $index);
	$index += scalar(keys%{ $currIon->{atoms} });
        print "Replacing $resname (residue# " . ($resId-$count) .") centered at " . PrintCoord($soluCoM) . " with " . $currIon->{NAME} . "\n";
    }

}

sub updateIonCoords {
    my ($ion, $soluCoM, $startIndex, $resId) = @_;
    my ($ionCoM, $i, $offset, $newIon, $j);

    $ionCoM = CoM($ion->{atoms});
    for $i ("XCOORD", "YCOORD", "ZCOORD") {
	$offset->{$i} = $soluCoM->{$i} - $ionCoM->{$i};
    }

    $newIon->{atoms} = dclone($ion->{atoms});
    $newIon->{bonds} = dclone($ion->{bonds});
    for $i (values %{ $newIon->{atoms} }) {
	$i->{RESNUM} = $resId;
	$i->{INDEX} += $startIndex;
	for $j ("XCOORD", "YCOORD", "ZCOORD") {
	    $i->{$j} += $offset->{$j};
	}
    }
    $newIon->{NAME} = "custom ion";
    return $newIon;
}

sub updateIonFields {
    my ($ionParms, $ionPos, $ionIndex, $ionResID) = @_;
    my (%ION, $newIon);
    
    %ION = %{ $ionPos };
    $ION{INDEX} = $ionIndex;
    $ION{ATMNAME} = $ionParms->{atoms}{1}{ATMNAME};
    $ION{RESNAME} = $ionParms->{atoms}{1}{ATMNAME};
    $ION{RESNUM} = $ionResID;
    $ION{FFTYPE} = $ionParms->{atoms}{1}{FFTYPE};
    $ION{NUMBONDS} = 0;
    $ION{LONEPAIRS} = 0;
    $ION{CHARGE} = $ionParms->{charge};
    $ION{RADII} = $ionParms->{atoms}{1}{RADII};
    $ION{LABEL} = "HETATM";
    $newIon->{atoms}{1} = \%ION;
    $newIon->{bonds}{1} = ();
    $newIon->{NAME} = $ionParms->{atoms}{1}{ATMNAME};
    return $newIon;
}

sub determineIonPlacement {
    my ($ions, $solvEnergies) = @_;
    my (@tmp, %CHARGE, $i);

    $ions->{1}{placed} = 0;
    $ions->{2}{placed} = 0;

    for $i (keys %{ $solvEnergies }) {
	$CHARGE{$solvEnergies->{$i}} = $i;
    }

    if ($ions->{1}{charge} > 0) {
	@tmp = sort numerically keys %CHARGE;
    } else {
	@tmp = reverse sort numerically keys %CHARGE;
    }

    while(@tmp) {
	last if ($ions->{1}{placed} >= $ions->{1}{total});
	$i = shift @tmp;
	push @{ $ions->{1}{list} }, $CHARGE{$i};
	$ions->{1}{placed}++;
    }
    if (keys %{ $ions->{2}{atoms}{1} }) {
	while (@tmp) {
	    last if ($ions->{2}{placed} >= $ions->{2}{total});
	    $i = pop @tmp;
	    push @{ $ions->{2}{list} }, $CHARGE{$i};
	    $ions->{2}{placed}++;
	}
    } else {
	delete $ions->{2};
    }
}

sub computeElecEnergyPerAtom {
    my ($atoms, $pqrFile) = @_;
    my ($apbsCmd, $isValid, %ENERGIES, $molID);
    
    $apbsCmd = "/home/yjn1818/codes/bin/coulomb -e ./_tmp.pqr";
    open COULCMD, "$apbsCmd |" or die "ERROR: Cannot run cmd $apbsCmd: $!\n";
    while (<COULCMD>) {
	chomp;
	if ($_ =~ /Atom\s+(\d+):\s+Energy\s+\=\s+(\-?\d+\.\d+E?.?\d*)/) {
	    if (exists($atoms->{$1}{IS_SOLVENT})) {
		$molID = ${ $atoms->{$1}{MOLECULEID} };
		$ENERGIES{$molID} += $2;
		$atoms->{$1}{ENERGY} = $2;
	    }
	}
    }
    die "ERROR: No valid energies obtained from $apbsCmd!\n" if (! %ENERGIES);

    system("rm -fr _tmp.pqr io.mc");
    return \%ENERGIES;
}

sub calcElecPot {
    my ($atoms, $bonds, $box, $headers) = @_;
    my ($pqrAtoms, $center, $sEng, $i);

    $pqrAtoms = dclone($atoms);
    $center = CoM($atoms);
    for $i ("X", "Y", "Z") {
	$center->{"${i}COORD"} = -($box->{$i}{len}/2 - $box->{$i}{lo} - $center->{"${i}COORD"});
    }
    &CenterSystem($pqrAtoms, $center, ());
    &createPQR($pqrAtoms, $bonds, "_tmp.pqr",$headers);
    $sEng = computeElecEnergyPerAtom($atoms, "_tmp.pqr");
    return $sEng;
}

sub getIonParms {
    my ($ionInfo, $parms) = @_;
    my ($counter, $element, $radii, $i, $count, $valid);

    if(exists($ionInfo->{1}{custom})) {
	delete $ionInfo->{2};
	return 0;
    }    

    for $counter (keys %{ $parms->{"ATOMTYPES"} }) {
        $element = $parms->{"ATOMTYPES"}{$counter}{"ATOM"};
        if ($element eq $ionInfo->{1}{atoms}{1}{"ATMNAME"}) {
            $ionInfo->{1}{atoms}{1}{"FFTYPE"} = $counter;
            &addRadii(\%{ $ionInfo->{1}{atoms}{1} }, $parms);
        }
        if (exists($ionInfo->{2}) and $element eq $ionInfo->{2}{atoms}{1}{"ATMNAME"}) {
            $ionInfo->{2}{atoms}{1}{"FFTYPE"} = $counter;
            &addRadii(\%{ $ionInfo->{2}{atoms}{1} }, $parms);
        }
    }

    if (! exists($ionInfo->{1}{atoms}{1}{"FFTYPE"})) {
        die "ERROR: The " . $ionInfo->{1}{atoms}{1}{"ATMNAME"} .
        " ion does not have any parameters in the provided forcefield\n";
    } elsif (exists($ionInfo->{2}) and ! exists($ionInfo->{2}{atoms}{1}{"FFTYPE"})) {
        die "ERROR: The " . $ionInfo->{2}{atoms}{1}{"ATMNAME"} .
        " ion does not have any parameters in the provided forcefield\n";
    }
    delete $ionInfo->{2} if (! keys %{ $ionInfo->{2} });
}

sub getIonOpts {
    my ($ionInfo, $BOX, $num_ions, $charge, $watCount) = @_;
    my ($mol_vol, $dim, $molar);

    $mol_vol = 1;
    for $dim ("X", "Y", "Z") {
        $mol_vol = $mol_vol * ($BOX->{$dim}{"hi"} - $BOX->{$dim}{"lo"});
    }
    printf "CELL VOLUME %.3G A^3", $mol_vol;
    if ($watCount) {
        $mol_vol = 18 * $watCount/.6023; # vol in cm^3
        printf "\nADJUSTED CELL MOLAR VOLUME based on 1g/cm3 WATER density %.3G cm^3", $mol_vol;
    }
    $mol_vol *= 6.0221415E-7;  #Angstroms^3 to meter^3
    $mol_vol *= 1000; #meter^3 to liter
    printf "(%.3G liter)\n", $mol_vol;
    if(! $ionInfo->{1}{custom}) {
	$ionInfo->{1}{charge} = $ionInfo->{1}{"total"} = 0;
	$ionInfo->{2}{charge} = $ionInfo->{2}{"total"} = 0;
    }

    if ((abs($charge - sprintf("%.0f", $charge))) > TOLERANCE) {
        print "WARNING: System has non-integer net charge of $charge!\n";
    }
    $charge = sprintf("%.0f", $charge);

    if (! $ionInfo->{1}{custom} and $ionInfo->{1}{atoms}{1}{ATMNAME} =~ /Na|K/i) { # +1 charge
        $ionInfo->{1}{charge} = "+1";
    } elsif (! $ionInfo->{1}{custom} and $ionInfo->{1}{atoms}{1}{ATMNAME} =~ /cl/i) {
        $ionInfo->{1}{charge} = -1;
    } elsif(! $ionInfo->{1}{custom}) {
        $ionInfo->{1}{charge} = "+2";
    }

    if (IsDecimal($num_ions) and $num_ions > 0) { #placing a concentration
        $ionInfo->{1}{"total"} = sprintf("%.0f",($mol_vol * $num_ions));
        $molar = $num_ions;
    } elsif ($num_ions > 0 and ! IsDecimal($num_ions)) {
        $ionInfo->{1}{"total"} = $num_ions;
        $molar = $ionInfo->{1}{"total"} / ($mol_vol);
    } elsif ($num_ions == 0 and ! IsDecimal($num_ions)) { # neutralize
        if ($ionInfo->{1}{charge} != 1) {
            if (($charge % 2) == 1) {
                print "WARNING: Charge is not even but divalent cation chosen! System will have net charge\n";
            }
	}
	die "ERROR: Ion and system have same charge but attempting to neutralize!\n" if(($charge/$ionInfo->{1}{charge})>0);
        $ionInfo->{1}{total} = sprintf("%.0f",-$charge/$ionInfo->{1}{charge});
        $molar = $ionInfo->{1}{total} / ($mol_vol);
    } else {
        die "ERROR: Invalid option for ion/salt type and concentration\n";
    }

    if (exists($ionInfo->{2}{atoms}{1}{ATMNAME})) {
        $ionInfo->{2}{charge} = -1;
        $ionInfo->{2}{total} = $ionInfo->{1}{total} * $ionInfo->{1}{charge};
        delete $ionInfo->{2} if ($ionInfo->{2}{total} < 0);
    }

    if (! exists($ionInfo->{2}{atoms}{1}{ATMNAME}) and ! exists($ionInfo->{1}{custom})) {
        printf "Placing " . $ionInfo->{1}{total} . " molecule(s) of " . $ionInfo->{1}{atoms}{1}{ATMNAME} .
            $ionInfo->{1}{charge} . " ion (%.3G molar)\n", $molar;
    } elsif (exists($ionInfo->{1}{custom})) {
	printf "Placing " . $ionInfo->{1}{total} . " molecule(s) of custom ion (" . $ionInfo->{1}{charge} . "e) " .
	    "(%.3G molar)\n", $molar;
    } else {
        printf "Placing " . $ionInfo->{1}{total} . " molecule(s) of " . $ionInfo->{1}{atoms}{1}{ATMNAME} .
            $ionInfo->{2}{atoms}{1}{ATMNAME} . " salt (%.3G molar)\n", $molar;
    }
    delete $ionInfo->{2} if (! exists($ionInfo->{2}{atoms}{1}{ATMNAME}));
}

sub updateSolvent {
    my ($atoms, $solvList) = @_;
    my ($i, $charge);

    $charge = 0;
    for $i (keys %{ $atoms }) {
	if (exists($solvList->{$i})) {
	    $atoms->{$i}{IS_SOLVENT} = 1;
	} else {
	    $charge += $atoms->{$i}{CHARGE};
	}
    }

    return $charge;
}

sub checkAtomTypes {
    my ($atoms, $parms) = @_;
    my ($i, $ffType);

    for $i (keys %{ $atoms }) {
	$ffType = $atoms->{$i}{FFTYPE};
	die "ERROR: Force field type $ffType not found in forcefield(s). Aborting\n"
	    if (! exists($parms->{ATOMTYPES}{$ffType}));
    	&addRadii($atoms->{$i}, $parms);
    }
}

sub addRadii {
    my ($atom, $parms) = @_;
    my ($i, $ffType);

    $ffType = $atom->{FFTYPE};
    $atom->{RADII} = $parms->{VDW}{$ffType}{$ffType}{1}{VALS}[1];
}

sub init {
    my (%OPTS, $forceFields, $ffType, $select, @tmp, $ION);
    getopt('bfwsin', \%OPTS);
    
    for ("i", "b", "f") {
	if (! exists($OPTS{$_})) {
	    &usage;
	    die "\n";
	}
    }
    print "Initializing...";
    ($bgfFile, $forceFields, $saveName, $select, $ionType, $ionConc) = 
	($OPTS{b}, $OPTS{f}, $OPTS{w}, $OPTS{s}, $OPTS{i}, $OPTS{n});
    FileTester($bgfFile);
    
    if (! defined($saveName)) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$//;
	$saveName .= "_ion.bgf";
    }

    ($ffFiles, $ffType) = ReadFFs($forceFields);

    $select = "NrWAT Sm3" if (! defined($select));
    @tmp = split /\s+/, $select;
    $solvSelect = GetSelections(\@tmp, 0);

    die "ERROR: Available ion types are - Na, Cl, Mg, K, Ca " .
        "and the corresponding salts (e.g. KCl)\n"
        if ($ionType !~ /Na|Cl|Mg|K|Ca/i and ! -e $ionType);

    $ionConc = 0 if (! defined($ionConc));
    if ($ionType !~ /^(Na|Mg|K|Ca|Cl)(Cl)?$/i and ! -e $ionType) {
        print "ERROR: Invalid ion/salt type $ionType\n";
        usage;
        die "\n";
    } elsif( ! -e $ionType) {
        $ION->{1}{atoms}{1}{"ATMNAME"} = $1;
        if ($2 and $1 ne $2) {
            $ION->{2}{atoms}{1}{"ATMNAME"} = $2;
        }
    } else {
	$ION->{1} = getIonType($ionType);
    }
    if ($ionConc eq "0" and length($ionType) > 2 and ! -e $ionType) {
        die "ERROR: Cannot use a salt to neutralize the system\n";
    }

    print "Done\n";

    return $ION;
}

sub getIonType {
    my ($ionFile) = $_[0];
    my ($ionAtoms, $ionBonds, $charge, $i, $data);

    ($ionAtoms, $ionBonds) = GetBGFFileInfo($ionFile,0,0);
    $charge = 0.0;
    for $i (keys %{ $ionAtoms }) {
	$charge += $ionAtoms->{$i}{CHARGE};
    }
 
    $charge += 0;
    die "ERROR: Charge read in from $ionFile ($charge) is not interger!\n" if($charge !~ /^\-?\d+$/);

    $data->{atoms} = $ionAtoms;
    $data->{bonds} = $ionBonds;
    $data->{charge} = $charge;
    $data->{total} = 0;
    $data->{custom} = 1;
    return $data;   
}
sub usage {
    print <<USAGE
usage: $0 -b bgf file -f force field(s) -i ion type -s [solvent selection] -n [\# ions] -w [save_name]
options:
    bgf_file: the name of the BGF formatted file
    force field: 1 (or more if in quotes) cerius2/mpsim/cmdf formatted forcefile for the bgf file
    ion_type: Na, Cl, Mg, K, Ca and the combination salts (e.g. NaCl) or add custom by giving bgf file
    solvent selection: Residues name (Nrxxx) or number of atoms per molecule (Imxxx) of the solvent
		       Defaults to NrWAT (name of residue = WAT) or Im3 (all molecules with 3 atoms)
		       for water
    [num_ions] : number of ions
                 0 - neutralize the system (default) - ions only
                 X - X number of ions/salt
                 X.x - X.x molar ion/salt
    [save_name]: name of output BGF file (optional)
USAGE
}

