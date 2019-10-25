#!/usr/bin/perl -w

use FindBin ();
use lib "$FindBin::Bin";

use strict;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use Math::GSL::Linalg::SVD;
use IO::Handle;

use Packages::General qw(FileTester TrjSelections CoM CrossProduct
			 DotProduct VecLen GetRotationMatrix);
use Packages::FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use Packages::ManipAtoms qw(UnwrapAtoms ScaleAtoms ImageAtoms SelectAtoms BuildAtomSelectionString
			GetAtmData GetMols RotateAbout GetAbituraryRotationMatrix);

sub init;
sub showUsage;
sub countMolsInCNT;
sub getDirectionVector;
sub saveSnapshots;
sub insideCylinder;
sub getCapPts;
sub alignAtoms;

my ($lammpsTrj, $bgfFile, $saveName, $prefix, $align, $offset, $use_radii);
my ($soluteAtoms, $solventAtoms, $MOLS, $SELECT);
my ($printAtoms, $LAMMPSOPTS, $pStr, $tMass);
my ($ATOMS, $BONDS, $SOLUTE, $SOLVENT);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
$SOLUTE = SelectAtoms($soluteAtoms, $ATOMS);
die "ERROR: No valid solute atoms found!\n" if (! $SOLUTE);
$tMass = 12.01*scalar(keys%{ $SOLUTE} );
$SOLVENT = GetMols($ATOMS, $BONDS, SelectAtoms($solventAtoms, $ATOMS));
die "ERROR: No valid solvent atoms found!\n" if (! $SOLVENT);
print "Done\n";
&GetLammpsByteOffset($SELECT, $lammpsTrj, scalar keys %{ $ATOMS });
&GetLammpsTrjType($SELECT, $lammpsTrj, "atom", \%{ $LAMMPSOPTS });
open OUTDATA, "> $saveName" or die "ERROR: Cannot write to $saveName: $!\n";
OUTDATA->autoflush(1);
print OUTDATA "#frame num_mols density radii height";
print OUTDATA " eff_density" if (defined($offset));
print OUTDATA " avail_density" if (defined($use_radii));
print OUTDATA "\n";
$pStr = "Parsing LAMMPS trajectory $lammpsTrj...";
ParseLAMMPSTrj($ATOMS, $lammpsTrj, $SELECT, "atom", \&countMolsInCNT, $pStr, \*OUTDATA);
close OUTDATA or die "ERROR: Cannot close $saveName: $!\n";

sub countMolsInCNT {
    my ($DATA, $bgfInfo, $fileHandle) = @_;
    my ($BOX, $totAtms, @tmp, $j, $i, $MOLECULE, $CENTER, $density_tot);
    my ($index, $sMol, $soluCenter, $u, $count, $v, $M, $density_eff, $density_avail);
    my ($radii_sq, $len_sq, $vfactor, $len, $pt1, $pt2, $eff_radii);

    $BOX = ConvertLammpsBox($DATA->{"BOX BOUNDS"});
    $totAtms = $DATA->{"NUMBER OF ATOMS"}[0];
    if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
	UnwrapAtoms($DATA->{ATOMS}, $BOX, $LAMMPSOPTS->{scaled});
    }
    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $j (@tmp) {
        $CENTER->{$j} = $BOX->{$j}{CENTER} = $BOX->{$j}{len}/2;
    }
    $SOLUTE = GetAtmData($DATA->{ATOMS}, $SOLUTE);
    $soluCenter = CoM($SOLUTE);
    for $j (@tmp) {
        $BOX->{$j}{hi} = $BOX->{$j}{len};
        $BOX->{$j}{lo} = 0;
        for $i (keys %{ $DATA->{ATOMS} }) {
    	    $DATA->{ATOMS}{$i}{$j} += ($BOX->{$j}{CENTER} - $soluCenter->{$j});
	}
    }
    $soluCenter = CoM($SOLUTE);
    #get direction vector
    ($u,$v,$M) = getDirectionVector($SOLUTE, $soluCenter);
    $radii_sq = $u->[2]/$tMass;
    $eff_radii = (sqrt($radii_sq)-$offset) if (defined($offset));
    $len_sq = ((($u->[0]+$u->[1])/2)*12/$tMass)-6*$radii_sq;
    $len = sqrt($len_sq);
    ($pt1, $pt2) = getCapPts($soluCenter, $v, $len/2);
    $vfactor = 0.6023*3.1492*$len;
    for $index (keys %{ $SOLVENT }) {
        $MOLECULE = GetAtmData($DATA->{ATOMS}, $SOLVENT->{$index}{MEMBERS});
	$CENTER = CoM($MOLECULE);
	ImageAtoms($MOLECULE, $CENTER, $BOX);
	$CENTER = CoM($MOLECULE);
	next if (! insideCylinder($pt1, $pt2, $len_sq, $radii_sq, $CENTER));
	$sMol->{$index} = $SOLVENT->{$index};
	for $j (keys %{ $CENTER }) {
	    $sMol->{$index}{$j} = $CENTER->{$j};
	}
    }
    $count = scalar keys %{ $sMol };
    $density_tot = $count*18.01/$vfactor/$radii_sq;
    $density_eff = $count*18.01/$vfactor/$eff_radii**2 if (defined($offset));
    $density_avail = $count*18.01/$vfactor/$use_radii**2 if (defined($use_radii));

    printf $fileHandle "%s %s %0.5f %0.5f %0.5f", $DATA->{TIMESTEP}[0],
		$count,$density_tot,sqrt($radii_sq),$len;
    printf $fileHandle " %.5f", $density_eff if (defined($offset));
    printf $fileHandle " %.5f", $density_avail if (defined($use_radii));
    print $fileHandle "\n";
    &alignAtoms($DATA->{ATOMS}, $soluCenter, $v) if ($align);
    &saveSnapshots($DATA->{ATOMS}, $BOX, $sMol, "${prefix}." . $DATA->{TIMESTEP}[0] . ".bgf") if ($prefix);
}

sub insideCylinder {
    my ($pt1, $pt2, $len_sq, $radii_sq, $testpt) = @_;
    my ($i, $d, $pd, $dot, $dist_sq);

    $dot = $dist_sq = 0.0;
    for $i ("XCOORD", "YCOORD", "ZCOORD") {
	$d->{$i} = $pt2->{$i} - $pt1->{$i};
	$pd->{$i} = $testpt->{$i} - $pt1->{$i}; # vector from cylinder to point
	$dot += $pd->{$i}*$d->{$i}; # dot product
	$dist_sq += $pd->{$i}*$pd->{$i}; #distance squared to cylinder axis
    }
    return 0 if($dot > $len_sq or $dot < 0); #dot product lies outside caps
    $dist_sq -= $dot*$dot/$len_sq;
    return 0 if ($dist_sq > $radii_sq);
    return 1;
}

sub getCapPts {
    my ($com, $dVec, $len) = @_;
    my ($pt1, $pt2, $i, $count);

    $count = 0;
    for $i ("XCOORD", "YCOORD", "ZCOORD") {
	$pt1->{$i} = $com->{$i} + $dVec->[$count]*$len;
	$pt2->{$i} = $com->{$i} - $dVec->[$count]*$len;
	$count++;
    }
    
    return ($pt1, $pt2);
}

sub getDirectionVector {
    my ($atoms, $com) = @_;
    my ($i, $j, $bounds, $mass, $pos, $inertia, $Matrix, $eigen);

    $mass = 12.011;
    $eigen = Math::GSL::Linalg::SVD->new( { verbose=>0 });
    for $i (values %{ $atoms }) {
	$pos->{XCOORD} = $pos->{YCOORD} = $pos->{ZCOORD} = $pos->{r} = 0.0;
	for $j ("XCOORD", "YCOORD", "ZCOORD") {
	    $pos->{$j} = $i->{$j} - $com->{$j};
	}
	$inertia->{Ixx} += $mass*($pos->{YCOORD}**2 + $pos->{ZCOORD}**2);
	$inertia->{Ixy} -= $mass*($pos->{XCOORD}*$pos->{YCOORD});
	$inertia->{Ixz} -= $mass*($pos->{XCOORD}*$pos->{ZCOORD});
        $inertia->{Iyy} += $mass*($pos->{XCOORD}**2 + $pos->{ZCOORD}**2);
        $inertia->{Iyz} -= $mass*($pos->{YCOORD}*$pos->{ZCOORD}); 
        $inertia->{Izz} += $mass*($pos->{XCOORD}**2 + $pos->{YCOORD}**2);
    }
    $Matrix = [[$inertia->{Ixx},$inertia->{Ixy},$inertia->{Ixz}],
		[$inertia->{Ixy},$inertia->{Iyy},$inertia->{Iyz}],
		[$inertia->{Ixz},$inertia->{Iyz},$inertia->{Izz}]];
    $eigen->load_data({ data=>$Matrix });
    $eigen->decompose( { algorithm=>q{jacobi} } );
    my ($eigenVals, $eigenVecs, $vMat, $oMat) = $eigen->results;
    return ($eigenVals, [$eigenVecs->[0][2],$eigenVecs->[1][2],$eigenVecs->[2][2]],$eigenVecs);
}

sub saveSnapshots {
    my ($atoms, $box, $sMol, $sname) = @_;
    my ($i, $j, $headers, @tmp);

    $box->{1}{DATA} = 90;
    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (2 .. 4) {
        $box->{$i}{DATA} = ($box->{$tmp[($i - 2)]}{hi} - $box->{$tmp[($i - 2)]}{lo});
    }

    $headers = createHeaders($box, $sname);

    for $i (keys %{ $ATOMS }) {
	for $j (@tmp) {
	    $ATOMS->{$i}{$j} = $atoms->{$i}{$j};
	}
	$ATOMS->{$i}{CHAIN} = "A" if (exists($SOLUTE->{$i}));
	$ATOMS->{$i}{CHAIN} = "X" if (!exists($SOLUTE->{$i}));
    }
    for $i (keys %{ $sMol }) {
	for $j (keys %{ $sMol->{$i}{MEMBERS} }) {
	    $ATOMS->{$j}{CHAIN} = "I";
	}
    }
    &addHeader($ATOMS, $headers);
    &createBGF($ATOMS, $BONDS, $sname);
}

sub alignAtoms {
    my ($atoms, $com, $dVec) = @_;
    my ($i, $j, $axis, $rotvec, $sine, $cosine, $angle, $vec, $rotM);

    $axis->{XCOORD} = $axis->{YCOORD} = 0;
    $axis->{ZCOORD} = 1;
    $vec->{XCOORD} = $dVec->[0];
    $vec->{YCOORD} = $dVec->[1];
    $vec->{ZCOORD} = $dVec->[2];
    #first image com to origin
    for $i (keys %{ $atoms }) {
	for $j ("XCOORD", "YCOORD", "ZCOORD") {
	    $atoms->{$i}{$j} -= $com->{$j};
	}
    }

    $rotvec = CrossProduct($vec,$axis);
    $cosine = DotProduct($vec,$axis);
    $sine = VecLen($rotvec);
    $angle = atan2($sine,$cosine);

    #rotate system
    $rotM = getAbituraryRotationMatrix($rotvec,$sine,$cosine);
    &rotateAbout($atoms, $rotM);
    #now place com back to box center
    for $i (keys %{ $atoms }) {
        for $j ("XCOORD", "YCOORD", "ZCOORD") {
            $atoms->{$i}{$j} += $com->{$j};
        }
    }
}

sub getAbituraryRotationMatrix {
    my ($v, $s, $c) = @_;
    my ($rM,$vx,$vy,$vz,$scale,$tmp);
    $tmp = VecLen($v);
    $scale = 1/$tmp;
    ($vx,$vy,$vz) = ($v->{XCOORD}*$scale,$v->{YCOORD}*$scale,$v->{ZCOORD}*$scale);
    $rM = [
	   [ 1+(1-$c)*($vx*$vx-1), (1-$c)*$vx*$vy-$vz*$s, (1-$c)*$vx*$vz+$vy*$s],
	   [(1-$c)*$vx*$vy-$vz*$s,  1+(1-$c)*($vy*$vy-1), (1-$c)*$vy*$vz-$vx*$s],
	   [(1-$c)*$vx*$vz-$vy*$s, (1-$c)*$vy*$vz+$vx*$s,  1+(1-$c)*($vz*$vz-1)]
	  ];
    return $rM;
}

sub rotateAbout {
    my ($atoms, $rotM) = @_;
    my ($orig, $i, $j, $count);

    for $i (keys %{ $atoms }) {
	for $j ("XCOORD","YCOORD","ZCOORD") {
	    $orig->{$i}{$j} = $atoms->{$i}{$j};
	}
    }

    for $i (keys %{ $atoms }) {
	$atoms->{$i}{XCOORD} =  $orig->{$i}{XCOORD}*$rotM->[0][0]+
				$orig->{$i}{YCOORD}*$rotM->[0][1]+
				$orig->{$i}{ZCOORD}*$rotM->[0][2];
	$atoms->{$i}{YCOORD} =  $orig->{$i}{XCOORD}*$rotM->[1][0]+
				$orig->{$i}{YCOORD}*$rotM->[1][1]+
				$orig->{$i}{ZCOORD}*$rotM->[1][2];
	$atoms->{$i}{ZCOORD} =  $orig->{$i}{XCOORD}*$rotM->[2][0]+
				$orig->{$i}{YCOORD}*$rotM->[2][1]+
				$orig->{$i}{ZCOORD}*$rotM->[2][2];
    }
}

sub init {
    my (%OPTS, $selection, $watAtoms, $cntAtoms, @tmp, $usage);

    getopt('lstbwnpaor',\%OPTS);
    $usage = &showUsage;

    ($lammpsTrj,$bgfFile,$selection,$saveName,$watAtoms,$cntAtoms,$prefix,$align,$offset,$use_radii) = 
	($OPTS{l},$OPTS{b},$OPTS{t},$OPTS{s},$OPTS{w},$OPTS{n},$OPTS{p},$OPTS{a},$OPTS{o},$OPTS{r});

    for ($lammpsTrj, $bgfFile) {
	die "$usage\n" if (! defined($_));
    }
    print "Initializing...";

    FileTester($lammpsTrj);
    FileTester($bgfFile);
    $selection = "*" if (! defined($selection));
    $SELECT = TrjSelections($selection);

   
    $watAtoms = "resname eq 'WAT'" if(! defined($watAtoms));
    $solventAtoms = BuildAtomSelectionString($watAtoms, 0);

    $cntAtoms = "resname ne 'WAT'" if (! defined($cntAtoms));
    $soluteAtoms = BuildAtomSelectionString($cntAtoms, 0);
    if(defined($prefix)) {
	$align = 0 if(! defined($align) or $align !~ /^(1|yes)$/i);
	$align = 1 if($align =~ /^(1|yes)$/i);
    }
    if (! defined($saveName)) {
	$saveName = $prefix if (defined($prefix));
	$saveName = basename($lammpsTrj) if (! defined($prefix));
	$saveName =~ s/\.\w+$//;
	$saveName .= ".cntWat.dat";
    }
    undef($offset) if (defined($offset) and ($offset !~ /\d+\.?\d*$/ or ! $offset));
    undef($use_radii) if (defined($use_radii) and ($use_radii !~ /\d+\.?\d*$/ or ! $use_radii));
    print "Done\n";
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -l lammps trj -t trajectory selection\n" . 
		  "-s [save name] -n [cnt atoms] -w [wat atoms] -p [save bgf prefix]\n" . 
		  "-o [radii offset=0] -r [available radii] -a [align to z axis]\n" .
	"options:\n" . 
	"-b bgf file: the location of the bgf file. This is required\n" . 
	"-l lammps trj: the location of the lammps trajectory file. This is required.\n" .
	"-t trajectory selection: The number of frames to use. Can be a single integer, or several integers in quotes\n" .
	"\t\tTo specify a range specify it as :Ita-b:c, which will take frames a - b, every c. Specify multiple ranges\n" .
	"\t\tor a combination of ranges and single frames by enclosing them in quotes. Specify \"*\" for all frames\n" .
	"-n [cnt atoms: resname ne 'WAT']: atom selection for cnt\n" .
	"-w [wat atoms: resname eq 'WAT']: atom selection for water\n" .
	"-s [save name]: name of output file\n" .
	"-p [save bgf prefix=none]: prefix for saving bgf files. will be {prefix}_tstep.bgf\n" .
	"-o [radii offset]: adjust the calculated tube radii by this amount to get effective density\n" .
	"-r [available radii]: calculate density assuming this tube radii as well\n" . 
	"-a [align to z axis=no]: output bgf files will be rotated to align cnt to z axis if yes\n";
    return $usage;
}
