#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use Packages::AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use Packages::General qw(TrjSelections FileTester GetBondLength CoM);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use Packages::FileFormats qw(GetBGFFileInfo GetMOL2FileInfo GetPDBFileInfo);
use Packages::ManipAtoms qw(SelectAtoms BuildAtomSelectionString SplitAtomsByMol GetMols);
use Getopt::Std qw(getopt);
use IO::Handle;

sub init;
sub calcNumAtoms;
sub updateBGF;

my ($readFunc, $dataFile, $aSELECT, $fileType, $tSELECT, $trjFile, $saveFile);
my ($ATOMS, $BONDS, $HEADERS, $calcMol, $MOLS, $cATOMS);
my ($field, $pStr, $LAMMPSOPTS, $getByteOffset, $getSnapshots, $trjType, $OUTFILE);

$|++;

&init;

print "Parsing $fileType file $dataFile...";
($ATOMS, $BONDS, $HEADERS) = $readFunc->($dataFile);
&GetMols($ATOMS, $BONDS);
delete $ATOMS->{HEADER} if (exists($ATOMS->{HEADER}));
print "Done\n";
$pStr = "Selecting Relevant Atoms...";
if (! defined($trjFile)) {
    print "$pStr";
    ($cATOMS, $BONDS) = SelectAtoms($aSELECT, $ATOMS);
    $MOLS = SplitAtomsByMol($ATOMS, $cATOMS) if ($calcMol);
    &calcNumAtoms($cATOMS, $MOLS);
} else {
    $field = scalar keys %{ $ATOMS };
    $getByteOffset->($tSELECT, $trjFile, $field);
    if ($trjType == 2) {
        &GetLammpsTrjType($tSELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
        $field = "coord";
    }
    $getSnapshots->($ATOMS, $trjFile, $tSELECT, $field, \&calcNumAtoms, $pStr, $OUTFILE);
    close $OUTFILE;
}

sub calcNumAtoms {
    my ($atoms, $box, $frameNum, $fileHandle) = @_;

    if (defined($trjFile) and $trjType == 2) { #LAMMPS
        $frameNum = $atoms->{TIMESTEP}[0];
        $box = ConvertLammpsBox($atoms->{"BOX BOUNDS"});
        %{ $box->{X} } = %{ $box->{XCOORD} };
        %{ $box->{Y} } = %{ $box->{YCOORD} };
        %{ $box->{Z} } = %{ $box->{ZCOORD} };
        $atoms = $atoms->{ATOMS};
        if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
            UnwrapAtoms($atoms,  $box, $LAMMPSOPTS->{scaled});
        }
    } elsif (defined($trjFile)) {
        $box = ConvertAmberBox(\%{ $box });
    }

    if (!defined($trjFile)) {
	print "Done\nFound " . scalar(keys %{ $atoms }) . " atoms\n";
	print "Found " . scalar(keys %{ $MOLS }) . " molecules\n" if ($calcMol);
    } else {
        &updateBGF($ATOMS, $atoms, "XCOORD", "YCOORD", "ZCOORD");
	($cATOMS, $BONDS) = SelectAtoms($aSELECT, $ATOMS);
	$MOLS = SplitAtomsByMol($ATOMS, $cATOMS) if ($calcMol);
	print $OUTFILE "$frameNum " . scalar(keys %{ $cATOMS });
        print $OUTFILE " " . scalar(keys %{ $MOLS }) if ($calcMol);
	print $OUTFILE "\n";
    }
}

sub updateBGF {
    my ($atoms, $newAtoms, @fields) = @_;
    my ($i, $j);

    for $i (keys %{ $atoms }) {
	for $j (@fields) {
	    $atoms->{$i}{$j} = $newAtoms->{$i}{$j};
	}
    }
}

sub init {
    my (%OPTS, $atomSelect, $trjSel, $list, $i);

    getopt('ftsmlwry',\%OPTS);
    die "usage: $0 -f data file -t (file type [bgf|pdb|mol2]) -l (traj file) -y (traj type=lammps(default) or amber)-r (traj range) -s [atom selection] -m (calc mols = no) -w (write traj data)\n" 
	if (! exists($OPTS{f}));

    print "Initializing...";
    ($dataFile, $atomSelect, $fileType, $calcMol, $trjFile, $trjType, $trjSel, $saveFile) = 
	($OPTS{f}, $OPTS{s}, $OPTS{t}, $OPTS{m}, $OPTS{l}, $OPTS{y}, $OPTS{r}, $OPTS{w});
    FileTester($dataFile);
    $atomSelect = "index > 0" if (! defined($atomSelect));
    if (! defined($fileType)) {
	if ($dataFile =~ /.*\.(\w+)$/) {
	    $fileType = $1;
	    if ($fileType !~ /^(bgf|mol2|pdb)$/i) {
		$fileType = "bgf";
	    }
	} else {
	    $fileType = "bgf";
	}
    }
    $fileType = uc($fileType);
    
    if ($fileType =~ /bgf|mol2|pdb/i) {
	$readFunc = eval ('\&Get' . $fileType . 'FileInfo');
    } else {
	die "ERROR: Expected bgf|mol2|pdb while parsing filetype. Got $fileType\n";
    }
    
    $aSELECT = BuildAtomSelectionString($atomSelect);
    $calcMol = 0 if (! defined($calcMol) or $calcMol !~ /(1|yes)/i);
    $calcMol = 1 if ($calcMol =~ /(1|yes)/i);

    if (defined($trjFile) and -e $trjFile and -r $trjFile and -T $trjFile) {
        if (! defined($trjType)) {
            if ($trjFile =~ /\.lammps/) {
                $trjType = "lammps";
            } else {
                $trjType = "amber";
            }
        }

        if (lc($trjType) ne "lammps") {
            $trjType = 1;
            $getSnapshots = \&ParseAmberTrj;
            $getByteOffset = \&GetAmberByteOffset;
        } else {
            $trjType = 2;
            $getSnapshots = \&ParseLAMMPSTrj;
            $getByteOffset = \&GetLammpsByteOffset;
        }
        if (! defined($trjSel)) {
            $trjSel = "*";
        }
        $list = TrjSelections($trjSel);
        for $i (keys %{ $list }) {
            $tSELECT->{$i} = $list->{$i};
        }
        die "ERROR: No valid frames selected with selection $trjSel!\n"
            if (! keys %{ $tSELECT } and $trjSel ne "*");
        if (! defined($saveFile)) {
           $saveFile = basename($trjFile);
           $saveFile =~ s/\.\w+$//;
           $saveFile .= ".count.dat";
        }
 	open $OUTFILE, "> $saveFile" or die "ERROR: Cannot write to $saveFile: $!\n";
	$OUTFILE->autoflush(1);
	print $OUTFILE "#tstep numAtoms";
        print $OUTFILE " numMols" if ($calcMol);
	print $OUTFILE "\n";
    }

    print "Done\n";
    
}
