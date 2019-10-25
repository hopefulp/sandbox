#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts/";
}

use strict;
use Packages::MolData;
use Getopt::Std qw(getopt);
use Packages::AMBER qw(AmberLib);
use File::Basename qw(basename);

sub init;
sub showUsage;
sub getCellOpts;
sub getSolvOpts;
sub embedMols;
sub findLigand;
sub typeLigand;
sub removeHydrogens;
sub runLeap;
sub createMembrane;
sub convertWAT;
sub createLAMMPSfiles;
sub getAtomRange;
sub execCmd;
sub getCellInfo;
sub addIons;
sub removeTmpFiles;

my ($sfile, $solute, $tmp, $tmp1);
my ($memOpts, $solvOpts, $ffields, $LIBS, $leapPrep, $LIGAND, $soluRange);

$|++;
&init;
print "Gettting data from solute $sfile->{type} file $sfile->{name}...";
$solute->read($sfile->{name}, $sfile->{type});
print "Done\nSearching for ligand...";
$LIGAND = &findLigand($solute, $LIBS, $sfile);
print "Done\n";
$leapPrep = &typeLigand($solute, $LIGAND, $sfile) if (defined($LIGAND));
print "Assigning AMBER atom types to system...";
$sfile->{pdb} = &removeHydrogens($sfile->{name}) if ($sfile->{type} ne "pdb");
$sfile->{pdb} = $sfile->{name} if ($sfile->{type} eq "pdb");
print "Done\n";
$tmp = &runLeap($sfile->{pdb},$leapPrep);
print "Embedding system into $memOpts->{size} $memOpts->{type} membrane...";
$tmp1 = &createMembrane($memOpts) if ($memOpts->{type} ne "popc" or $memOpts->{size} ne "80 80");;
$sfile->{save} =~ s/\.\w+$/\.bgf/;
&embedMols($tmp, $tmp1, $ffields, $sfile->{save});
print "Done\nNeutralizing system...";
&addIons($sfile->{save}, $ffields);
print "Creating $sfile->{save}...";
&convertWAT($sfile->{save}, $solvOpts) if (defined($solvOpts));
print "Done\nCreating LAMMPS input and data files...";
&createLAMMPSfiles($sfile->{save}, $ffields);
print "Done\n";
&removeTmpFiles($sfile);

sub removeTmpFiles {
    my ($sfile) = $_[0];
    my ($prefix);

    print "Removing temporary files...";
    $prefix = $sfile->{name};
    $prefix =~ s/\.\w+$//;
    
    &execCmd("rm -fr ???.frcmod ???.prepi _tmplig.* ${prefix}_noH.* leaprc membrane.bgf _gpcrmembranemd.* leap.log _gpcrmembranecharge.dat ${prefix}_amber.bgf");
    print "Done\n";
}

sub addIons {
    my ($fileName, $ffields) = @_;
    my ($chargeCmd, $neutCmd, $charge, $rCharge, $delCharge, $addChargeCmd);

    $chargeCmd = "/home/yjn1818/scripts/bgfcharge.pl $fileName";
    $addChargeCmd = "/home/yjn1818/scripts/modifyAtomData.pl -s $fileName -w $fileName -a Ia1 -f \"CHARGE:";
    $neutCmd = "/home/yjn1818/scripts/addIons.pl -b $fileName -n 0 -f \"$ffields\" -w $fileName -i ";
    &execCmd($chargeCmd, "_gpcrmembranecharge.dat");
    open CHARGECMD, "_gpcrmembranecharge.dat" or die "ERROR: Cannot execute $chargeCmd: $!\n";
    while (<CHARGECMD>) {
	chomp;
	if ($_ =~ /Total Charge: (\-?\d+\.?\d*)/) {
	    $charge = $1;
	}
    }
    close CHARGECMD;
    die "ERROR: No valid data read while executing cmd $chargeCmd!\n" if (! defined($charge));
    $rCharge = sprintf("%.0f",$charge);
    $delCharge = $rCharge - $charge;
    if (abs($delCharge) > 0.00001) {
	$addChargeCmd .= "+" if ($delCharge > 0);
	$addChargeCmd .= $delCharge;
	&execCmd("${addChargeCmd}\"");
    }
    if ($rCharge > 0) {
	$neutCmd .= "Cl";
    } else {
	$neutCmd .= "Na";
    }
    &execCmd($neutCmd);
}

sub createLAMMPSfiles {
    my ($saveName, $ffields) = @_;
    my ($lammpsStr, $prefix);

    $prefix = basename($saveName);
    $prefix =~ s/\.\w+$//;
    
    $lammpsStr = "/home/yjn1818/scripts/createLammpsInput.pl -f \"$ffields\" -b $saveName -t gpcr -s $prefix";
    &execCmd($lammpsStr);
}

sub convertWAT {
    my ($filename, $opts) = @_;
    my ($convertStr, $i, $tmp, $j);

    print "Converting water model...";
    $convertStr = "/home/yjn1818/scripts/modifyAtomData.pl -s $filename -w $filename";
    for $i (keys %{ $opts }) {
	$tmp = "-a Ta${i} -f \"";
	for $j (keys %{ $opts->{$i} }) {
	    next if (! exists($opts->{$i}{$j}));
	    $tmp .= "${j}:$opts->{$i}{$j} ";
	}
	&execCmd("${convertStr} ${tmp} \"");
    }
}

sub createMembrane {
    my ($options) = $_[0];
    my ($CELL, $i, $DIMS, $replicateStr, $memFile, $shouldReplicate, $trimStr, $multiple);

    $memFile = "/home/yjn1818/Research/membranes/" . $options->{type} . "_80x80_amber_qeq_tip4p_solv.bgf";
    $replicateStr = "/home/yjn1818/scripts/replicate.pl -b $memFile -s ./membrane.bgf -d \"";
    $trimStr = "/home/yjn1818/scripts/trimCell.pl -b ./membrane.bgf -s _gpcrmembranemd.bgf -m 1 ";    
    $CELL = getCellInfo($options, $memFile);

    $shouldReplicate = 0;
    for $i ("x", "y") {
	$DIMS->{$i} = 1;
	$multiple = int($options->{$i}/$CELL->{$i}) + 1;
	$DIMS->{$i} = $multiple if ($multiple > 1);
	$shouldReplicate = 1 if ($multiple > 1);
	$replicateStr .= "$DIMS->{$i} ";
    }    
    $replicateStr .= "1\"";
    &execCmd($replicateStr) if ($shouldReplicate);
    &execCmd("cp $memFile ./membrane.bgf") if (! $shouldReplicate);
    $trimStr .= "-c \"$options->{x} $options->{y} $CELL->{z}\"";
    &execCmd($trimStr);
    return "_gpcrmembranemd.bgf";
}
sub runLeap {
    my ($pdbFile, $leapAdd) = @_;
    my ($leaprcLoc, $leapCmd, $molname, $convertCmd, $saveName);
    
    $molname = $pdbFile;
    $molname =~ s/\.\w+//;
    $saveName = $molname;
    $saveName =~ s/_noH//;
    $saveName .= "_amber.bgf";

    $leaprcLoc = "/home/yjn1818/amber/leaprc";
    $leapCmd = "/home/yjn1818/programs/ambertools/exe/tleap";
    $convertCmd = "/home/yjn1818/scripts/amber2bgf.pl ${molname}.prmtop ${molname}.inpcrd $saveName";
    &execCmd("cp $leaprcLoc ./");
    open MYLEAPRC, ">> leaprc" or die "Cannot write to leaprc: $!\n";
    print MYLEAPRC "gaff = loadamberparams gaff.dat\n";
    print MYLEAPRC $leapAdd if (defined($leapAdd));
    print MYLEAPRC "prot = loadpdb $pdbFile\n";
    print MYLEAPRC "saveamberparm prot ${molname}.prmtop ${molname}.inpcrd\n";
    print MYLEAPRC "quit\n";
    close MYLEAPRC;
    &execCmd("$leapCmd");
    print "Converting AMBER files to bgf...";
    &execCmd($convertCmd);
    print "Done\n";
    return $saveName;
}

sub removeHydrogens {
    my ($fileName) = $_[0];
    my ($cmd, $saveName); 
    $cmd = "/project/Biogroup/scripts/perl/bgf2pdb_noH.pl";
    $saveName = basename($fileName);
    $saveName =~ s/\.\w+$//;
    $saveName .= "_noH.pdb";

    print "Removing all hydrogens and converting to pdb...";
    &execCmd("${cmd} $fileName", $saveName);
    return $saveName;
}

sub typeLigand {
    my ($solu, $ligands, $filedata) = @_;
    my ($antechamberCmd, $i, $extractCmd, $resname, $atomRange, $convertCmd);
    my ($charge, $RES, $topo, $leapAdd, $saveName, $parmchkCmd, $resParm);

    print "Assigning atom types for ligands...";
    $antechamberCmd = "/exec/amber9/exe//antechamber -fo prepi -c am1 -at gaff -pf y";
    $extractCmd = "/home/yjn1818/scripts/getBGFAtoms.pl -b $filedata->{name} -s _tmplig.bgf -o";
    $convertCmd = "/home/yjn1818/scripts/bgf2mol2.pl -v 0 -b _tmplig.bgf -a 0";
    $parmchkCmd = "/exec/amber9/exe//parmchk -f prepi";
    for $resname (keys %{ $ligands }) {
	$charge = 0;
	$RES = uc $resname;
	$RES = substr($resname, 0, 3) if (length($resname) > 3); 
        $resParm = lc $RES . "Parm";
	$topo = "${resname}.res";
	$saveName = "${resname}.prepi";
	$atomRange = getAtomRange($ligands->{$resname});
	for $i (values %{ $ligands->{$resname} }) {
	    $charge += $i->charge;
	}
	$charge = 0 if (abs($charge) < 0.01);
	&execCmd("${extractCmd} \"$atomRange\"");
	&execCmd($convertCmd);
	&execCmd("${antechamberCmd} -fi mol2 -i _tmplig.mol2 -nc $charge -rn $RES -rf $topo -o $saveName");
        &execCmd("${parmchkCmd} -i ${RES}.prepi -o ${RES}.frcmod");
	$leapAdd .= "loadamberprep $saveName\n";
        $leapAdd .= "$resParm = loadamberparams ${RES}.frcmod\n";
    }
    
    print "Done\n";
    return $leapAdd;
}


sub findLigand {
    my ($soluFile, $amberLibs) = @_;
    my ($resname, $LIGANDS, $count, $i, @atomList, $atomId);

    $count = 0;
    for $i (keys %{ $soluFile->shash->{resid} }) {
        @atomList = keys %{ $soluFile->shash->{"resid"}{$i} };
        next if (! @atomList);
        $atomId = shift @atomList;
        next if (! $atomId or ! $soluFile->shash->{"resid"}{$i}{$atomId});
        $resname = $soluFile->shash->{"resid"}{$i}{$atomId}->resname;	
        $resname = "HIE" if ($resname =~ /his|hse|hdd/i);
        $resname = "CYS" if ($resname =~ /cyx/i);
	if (! exists($amberLibs->{uc $resname}) and ! exists($LIGANDS->{$resname})) {
	    $LIGANDS->{$resname} = $soluFile->shash->{resid}{$i};
	    $count++;
	}
    }
    print "none found..." if (! $count);
    print "found $count..." if ($count);

    return $LIGANDS;
}

sub embedMols {
    my ($soluFile, $solvFile, $ffields, $saveName) = @_;
    my ($embedStr);

    $embedStr = "/home/yjn1818/scripts/embedMolecule.pl -m $solvFile -s $soluFile -f \"$ffields\" -c com -w $saveName";
    &execCmd($embedStr);
    &execCmd("rm -fr _out _solv_1cell.bgf _solv_replicate.bgf _solv_trim.bgf");
}

sub createSolventBox {
    my ($solv, $solu, $cell) = @_;
    my ($i, $replicate, $blen, $trim, $box, $bmin, $smin, $offset, $map, $solvBox, $cellScale);

    $map = ({ "a" => "x", "b" => "y", "c" => "z"});
    if (defined($cell->{density})) { #compress/expand the solvent cell to the new density. assume 1 g/cm3
	$cellScale = 1/($cell->{density}**(1/3)); #
	for $i ("x", "y", "z") {
	    $solv->stressCell($i, $cellScale);
	}
    }
    $smin = $solv->getExtrema("min");
    %{ $solvBox } = %{ $solv->vdwbox };
    %{ $solvBox } = %{ $solv->cell } if ($solv->cell->{valid});
    %{ $box } = %{ $solu->vdwbox };
    $replicate = "/home/yjn1818/scripts/replicate.pl -b ./_solv_1cell.bgf -d \"";
    for $i ("a", "b", "c") {
	$box->{$i}{max} += $cell->{cell}{$i}{max} 
		if (exists($cell->{cell}) and exists($cell->{cell}{$i}) and exists($cell->{cell}{$i}{max}));
        $box->{$i}{min} -= $cell->{cell}{$i}{min} 
		if (exists($cell->{cell}) and exists($cell->{cell}{$i}) and exists($cell->{cell}{$i}{min}));
	$box->{$i}{len} = $box->{$i}{max} - $box->{$i}{min};
	$bmin->{$i} = $box->{$i}{min};
	$blen .= "$box->{$i}{len} ";
	$replicate .= sprintf("%.0f ", (($box->{$i}{len}/$solvBox->{$i}))+1);
    }
    $replicate .= "\" -s _solv_replicate.bgf";
    #move the solvent box to the solute minima
    for $i ("a", "b", "c") {
	$offset->{ $map->{$i }} = $bmin->{$i} - $smin->{$i};
    }
    $solv->moveMol("all", $offset);
    $solv->write("_solv_1cell.bgf", "bgf");
    #replicate the cell by the replication vector calculated above
    die "ERROR while executing \"$replicate\"\n" if (system("${replicate} >& _out.dat"));
    # remove all molecules outside the solute (inflated) cell
    $trim = "/home/yjn1818/scripts/trimCell.pl -b _solv_replicate.bgf -c \"$blen\" -s _solv_trim.bgf -m 1 -o 2";
    die "ERROR while executing \"$trim\"\n" if (system("${trim} >& _out.dat"));
    undef($solv);
}

sub init {
    my (%OPTS, $ffStr, $memTypeStr, $solvTypeStr, $periodStr, $memSizeStr);

    getopt('idmsw', \%OPTS);
    die &showUsage . "\n" if (! exists($OPTS{i}));

    print "Initialzing...";
    ($sfile->{name}, $memSizeStr, $memTypeStr, $solvTypeStr, $sfile->{save}) = 
	($OPTS{i}, $OPTS{d}, $OPTS{m}, $OPTS{s}, $OPTS{w});
    $solute =  Packages::MolData->new();
    $solute->testFile($sfile->{name});
    if ($sfile->{name} =~ /\.(\w+)$/) {
	$sfile->{type} = lc $1;
    } else {
	die "ERROR: Cannot determine file type from $sfile->{name}!\n";
    }
    $memTypeStr = "popc" if (! defined($memTypeStr) or $memTypeStr !~ /^(popc|pope)$/i);
    $memSizeStr = "80 80" if (! defined($memSizeStr));
    $memOpts = getCellOpts($memTypeStr, $memSizeStr);
    $ffields = "/home/yjn1818/ff/GAFF.ff /home/yjn1818/ff/AMBER03.ff /home/tpascal/ff/WAT/tip4ew.ff";
    ($solvOpts, $ffields) = getSolvOpts($solvTypeStr) if (defined($solvTypeStr));
    $sfile->{save} = $solute->getFileName($sfile->{name}) if (! defined($sfile->{save}));
    $LIBS = &AmberLib;
    system("rm -fr gpcrlammpsmd.log");
    $ENV{AMBERHOME} = "/home/yjn1818/programs/ambertools";
    print "Done\n";
}

sub getCellOpts {
    my ($memType, $memStr) = @_;
    my ($CELL);
    $CELL->{type} = $memType;
    if ($memStr =~ /^\s*(\d+\.?\d*)\s+(\d+\.?\d*)/) {
	($CELL->{x}, $CELL->{y}) = ($1, $2);
    } else {
	print "invalid membrane size \"$memStr\" encountered.. defaulting to 80 x 80...";
	($CELL->{x}, $CELL->{y}) = (80, 80);
    }
    $CELL->{size} = "$CELL->{x}x$CELL->{y}";
    return $CELL;
}

sub getSolvOpts {
    my ($solventStr) = lc $_[0];
    my ($SOLVENT, $ffields);
    
    $ffields = "/home/yjn1818/ff/GAFF.ff /home/yjn1818/ff/AMBER03.ff ";
    
    if ($solventStr =~ /tip3/) {
	$SOLVENT->{OW}{CHARGE} = "+0.2144";
	$SOLVENT->{HW}{CHARGE} = "-0.1572";
	if ($solventStr =~ /charmm/) {
	    $ffields .= "/home/yjn1818/ff/WAT/tip3_charmm.ff";
	} else {
	    $ffields .= "/home/yjn1818/ff/WAT/tip3ew.ff";
	}
    } elsif ($solventStr =~ /spc/) {
	$SOLVENT->{OW}{CHARGE} = "+0.2008";
	$SOLVENT->{HW}{CHARGE} = "-0.1004";
	$ffields .= "/home/yjn1818/ff/WAT/spcew.ff";
    } elsif ($solventStr =~ /f3c/) {
	$SOLVENT->{OW}{CHARGE} = "+0.2284";
	$SOLVENT->{HW}{CHARGE} = "-0.1142";
	$SOLVENT->{OW}{FFTYPE} = "O_3F";
	$SOLVENT->{HW}{FFTYPE} = "H_F";
	$ffields .= "/home/yjn1818/ff/WAT/f3cew.ff";
    }
    return ($SOLVENT, $ffields);
}

sub getAtomRange {
    my ($list) = $_[0];
    my (@atomids, $i, $prev, $start, $range);

    @atomids = sort {$a<=>$b} keys %{ $list };
    $i = 1;
    $prev = $start = $atomids[0];;
    while ($i <= $#atomids) {
	if (($atomids[$i] - $prev) > 1) {
	    $range .= ":Ia${start}-${prev} ";
	    $start = $atomids[$i];
	}
	$prev = $atomids[$i];
	$i++;
    }
    $range .= ":Ia${start}-${prev} ";
    return $range;
}

sub getCellInfo {
    my ($options, $memFile) = @_;
    my ($crystxCmd, $CELL);

    #$crystxCmd = "grep crystx $memFile";
    #&execCmd($crystxCmd, "_crystx.dat");
    open CRYSTXCMD, "grep -i crystx $memFile |" or die "ERROR while executing $crystxCmd |: $!\n";
    while (<CRYSTXCMD>) {
	chomp;
	if ($_ = ~/(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)/) {
	    ($CELL->{x}, $CELL->{y}, $CELL->{z}) = ($1, $2, $3);
	}
    }
    close CRYSTXCMD;
    die "ERROR: No valid cell info found while search $memFile!\n" if (! defined($CELL));
    return $CELL
}

sub execCmd {
    my ($cmdStr, $output) = @_;

    if (! defined($output)) {
	$cmdStr .= ">> gpcrlammpsmd.log";
	$output = "gpcrlammpsmd.log";
    } else {
	$cmdStr .= ">& $output";
    }
    system("echo \"\ncommand:'$_[0]'\" >> $output");
    if (system($cmdStr)) {
        die "ERROR: Check grcrlammpsmd.log\n";
    }
}

sub showUsage {
    my $usage = <<DATA;
usage: $0 -i input structure -m (membrane type) -d (membrane xy size) -s (solvent type) -w (savename)
OPTIONS:
	-i input structure:  REQUIRED. The following file formats are supported: BGF, PDB, MSI, MOL2
	-m membrane type:    OPTIONAL. Either POPC or POPE. Defaults to POPC.
	-d membrane xy size: OPTIONAL. The size of the membrane to construct. Needs to be at least as large as the protein
	                     Will default to 80 x 80.
	-s solvent type:     OPTIONAL. The following predetermined (equilibrated) solvents are available. Default F3C.
			     TIP3: the original Jorgenson TIP3 water model (rigid hydrogens)
			     TIP4: TIP4P with massless pseudo-atom
			     TIP3_CHARMM: TIP3 water model as implemented in CHARMM
			     F3C: F3C water model (no rigid hydrogens)
			     SPC: SPC water model
	-w savename:         OPTIONAL. Will assume \$prefix_mod.\$suffix if not specified.
DATA

    return $usage;
}
