#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
no warnings;
use Packages::General qw(FileTester LoadElements Permutate GetSoluteAtoms CRadDegrees IsDecimal LoadFFs);
use Packages::FileFormats qw(GetBGFFileInfo GetBondList);
use Packages::CERIUS2 qw(ReadFFs);
use Packages::BOX qw(GetBox);
use Packages::ManipAtoms qw(FindElement GetMols);
use Packages::LAMMPS qw(CreateInputFile);
use Getopt::Std qw(getopt);

sub init;
sub usage;
sub checkAtomTypes;
sub updateTorsionList;
sub getValenceParms;
sub numerically { ($a<=>$b) }
sub getLeastX;
sub findValenceType;
sub getPermutations;
sub getValParent;
sub updateParmIndex;
sub createDatFileHeader;
sub printValence;
sub sortParmByIndex;
sub getTorsionWeightFactor;
sub removeEmptyParms;
sub is5MemberRing;
sub getAtmList;
sub getBonds;
sub getCorrectInversionIndices;
sub getTorsionScalingFactor;
sub getNewParm;
sub determineIfHybrid;
sub findVDWEntry;
sub getPairMix;
sub getImage;
sub createLammpsClusterScript;
sub findDuplicateParms;
sub getElement;
sub getTIP4Popts;
sub getInputTypes;
sub setOpts;

my ($bgfFile, $FF, $suffix, %OPTS, $inputType, $ScaleTorsionList, $reaxFF, $opts);
my ($ATOMS, $CONS, $HEADERS, $FILES, $PARMS, %ERRORS, $i, $BOX, $ffType);
my ($BONDS, $TORSIONS, $INVERSIONS, $ANGLES, @atmIndices, $totAtms);
my ($writeLabels, $MOLS);

$|++;
my ($start) = time();
$FILES = &init;
$PARMS = LoadFFs($FILES, 2);
print "Step 3: Parsing BGF file $bgfFile...";
($ATOMS, $CONS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$MOLS = GetMols($ATOMS, $CONS);
$BOX = GetBox($ATOMS, $PARMS, $HEADERS);
&checkAtomTypes($ATOMS, $PARMS->{ATOMTYPES}, $ffType);
$CONS = () if ($ffType == 5); # remove all valence interactions for reax/3 body force fields
&findDuplicateParms($PARMS);
&updateTorsionList($PARMS->{TORSIONS});
&setOpts($ATOMS, $PARMS, $opts);
print "Done\n";
&printErrors(\%ERRORS, 1) if (keys %ERRORS);
print "Step 4: Determining valence list from connectivities...";
@atmIndices = sort numerically keys %{ $ATOMS };
&getValenceParms($ATOMS, $CONS, $PARMS, \@atmIndices);
&updateParmIndex($PARMS, $ffType);
&getTorsionScalingFactor($TORSIONS);
print "Done\n";
open DATFILE, "> data.${suffix}" || die "ERROR: Cannot create data.${suffix}: $!\n";
&createDatFileHeader($PARMS, $BOX, \*DATFILE);
&printAtoms($ATOMS, $BOX, \*DATFILE);
for $i ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
    &printValence(eval('$' . $i), \*DATFILE, $i);
}
close DATFILE;
print "Step 5: Writing data file data.${suffix}...";
print "will use single inversion..." if ($PARMS->{PARMS}{single_inversion});
&addLammpsParms($PARMS, $ATOMS, $MOLS, $ANGLES);
&getTIP4Popts($PARMS, $ATOMS, $BONDS) if (exists($PARMS->{PARMS}{is_tip4p}));
print "Done\nStep 6: Creating $PARMS->{PARMS}{INPUTTYPE} input files in.${suffix}...";
&CreateInputFile($PARMS);
print "Created in.${suffix}_singlepoint...Done\n";
print "Step 7: Creating LAMMPS cluster script file ${suffix}.lammps.pbs...";
&createLammpsClusterScript($FILES);
print "Done\n";
my ($end) = time();
printf "Elapsed time: %.3f secs\n", ($end - $start);
&printErrors(\%ERRORS, 0) if (keys %ERRORS);

sub printValence {
    my ($valence, $datFile, $header) = @_;
    my ($i, $j, $count);
    
    $header = "dihedrals" if ($header eq "TORSIONS");
    $header = "impropers" if ($header eq "INVERSIONS");
    
    $header = uc(substr($header, 0, 1)) . lc (substr($header,1,length($header)));
    print "Step 5: Writing data file...$header\r";
    print $datFile "\n$header\n";
    for $i (1 .. $valence->{counter}) {
	printf $datFile "\n%8d %8d ", $i, $valence->{LIST}{$i}{DATA}{INDEX};
	$count = 0;
	for $j (@{ $valence->{LIST}{$i}{ATOMS} }) {
	    printf $datFile "%8d ", $j;
	    $count++;
	}
    }
    print $datFile "\n";
}
	
sub createDatFileHeader {
    my ($parms, $box, $datFile) = @_;
    my ($index, $curr, $parm, $DATA, $count, $j, $valence);
    my ($a, $b, $c, $cosA, $cosB, $cosC, $parmHybrid);

    $parms->{PARMS}{isTriclinic} = 0;
    print "Step 5: Writing data file...header\r";
    print $datFile "Created by $0 on " . scalar(localtime) . "\n\n";
    printf $datFile "%12d atoms\n", $totAtms;
    for $index ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
	$curr = eval('$' . $index);
	$parm = lc($index);
	$parm = "dihedrals" if ($parm eq "torsions");
	$parm = "impropers" if ($parm eq "inversions");

	printf $datFile "%12d $parm\n", $curr->{counter}
    }
    printf $datFile "\n%12d atom types\n", $parms->{ATOMTYPES}{counter};
    for $index ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
	$parm = lc(substr($index,0,-1));
	$parm = "dihedral" if ($parm eq "torsion");
	$parm = "improper" if ($parm eq "inversion");
	printf $datFile "%12d $parm types\n", $parms->{$index}{counter};
    }
    print "Step 5: Writing data file...box\r";
    for $index ("X","Y","Z") {
	printf $datFile "\n         %10.6f %10.6f ", $box->{$index}{"lo"}, $box->{$index}{"hi"};
	print $datFile lc($index) . "lo " . lc($index) . "hi";
	$parms->{PARMS}{isTriclinic} = 1 if ($box->{$index}{angle} != 90);
    }
    
    if ($parms->{PARMS}{isTriclinic}) { #write xy xz yz
	$a = $box->{X}{len};
	$b = $box->{Y}{len};
	$c = $box->{Z}{len};
	$cosA = cos(CRadDegrees($box->{X}{angle}, 0));
	$cosB = cos(CRadDegrees($box->{Y}{angle}, 0));
	$cosC = cos(CRadDegrees($box->{Z}{angle}, 0));

	printf $datFile "\n         %10.6f %10.6f %10.6f xy xz yz",
        $b*$cosA,
        $c*$cosC,
	$c*($cosB-($cosA*$cosC))/sqrt(1-($cosA*$cosA));
    }
    print "Step 5: Writing data file...coeffs\r";
    ($DATA, $count) = sortParmByIndex($parms->{VDW}, 1);
    print $datFile "\n\nMasses\n\n";
    for $index (1 .. $count) {
	printf $datFile "%5d %8.4f", $index, $DATA->{$index}{DATA}{MASS};
	printf $datFile "%-10s",' # ' . $DATA->{$index}{DATA}{NAME} if ($writeLabels);
	print $datFile "\n";
    }
    ($DATA, $count) = sortParmByIndex($parms->{VDW}, 1);
    if ($ffType != 5 && (scalar(keys %{ $parms->{VDW}{TYPE} }) == 1 && $DATA->{1}{Lammps}{name} ne "sw")) {
	print $datFile "\nPair Coeffs\n";
    
	for $index (1 .. $count) {
	    printf $datFile "\n%5d ", $index;
	    if (scalar(keys %{ $parms->{VDW}{TYPE} }) > 1) {
		printf $datFile "%15s ", $DATA->{$index}{Lammps}{name};
	    }
	    for $j (@{ $DATA->{$index}{VALS} }) {
		last if (! defined($j) or $j eq "" or $j =~ /\#/);
		if (IsDecimal($j)) {
		    printf $datFile "%15.8f ", $j;
		} else {
		    printf $datFile "%15d ", $j;
		}
	    }
	}
    }

    $valence = "";   
    for $curr ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
	next if ($parms->{$curr}{counter} == 0);
	$parm = uc(substr($curr,0,1)) . lc(substr($curr,1,-1));
	$parm = "Dihedral" if ($parm eq "Torsion");
	$parm = "Improper" if ($parm eq "Inversion");
	($DATA, $count) = sortParmByIndex($parms->{$curr});
		printf $datFile "%15s ", $DATA->{$index}{Lammps}{name};
	if (scalar(keys %{ $parms->{$curr}{TYPE} }) == 1) {
	    print $datFile "\n\n${parm} Coeffs\n";
	    for $index (1 .. $parms->{$curr}{counter}) {
		printf $datFile "\n%5d ", $index;
		if (! $parms->{PARM}{single_inversion} and $curr eq "INVERSIONS" and 
		    $DATA->{$index}{Lammps}{name} eq "umbrella") {
		    #$DATA->{$index}{VALS}[0] /= 3;
		}
		for $j (@{ $DATA->{$index}{VALS} }) {
		    #last if (! defined($j) or $j eq "" or $j =~ /\#/);
		    printf $datFile "%15g ", $j;
	        }
	        print $datFile '# ' . $DATA->{$index}{KEY} if ($writeLabels);
	    }
	} else {
	    $valence = lc $curr;
	    chop $valence;
	    for $index (1 .. $parms->{$curr}{counter}) {
		$parms->{HYBRID_VALENCE} .= sprintf("${valence}_coeff\t%5d %-15s",$index,$DATA->{$index}{Lammps}{name});
		if (! $parms->{PARM}{single_inversion} and $curr eq "INVERSIONS" and
                    $DATA->{$index}{Lammps}{name} eq "umbrella") {
                    $DATA->{$index}{VALS}[0] /= 3;
                }
                for $j (@{ $DATA->{$index}{VALS} }) {
                    #last if (! defined($j) or $j eq "" or $j =~ /\#/);
                    if (IsDecimal($j)) {
                        $parms->{HYBRID_VALENCE} .= sprintf("%15.8f ", $j);
                    } elsif($j =~ /[a-zA-Z]/) {
                        $parms->{HYBRID_VALENCE} .= sprintf("%15s ", $j);
                    } else {
                        $parms->{HYBRID_VALENCE} .= sprintf("%12d ", $j);
                    }
                }
                $parms->{HYBRID_VALENCE} .= sprintf('# ' . $DATA->{$index}{KEY}) if ($writeLabels);
		$parms->{HYBRID_VALENCE} .= "\n";
            }
	    $parms->{HYBRID_VALENCE} .= "\n";
        }
    }
    printf $datFile "\n\n";
}

sub printAtoms {
    my ($atoms, $box, $datFile) = @_;
    my ($counter, $type_id, $atm_name, $fmt, $out_string, $index, $dim, %IMAGE);

    $fmt = "%8d %8d %8d %11.5f %10.5f %10.5f %10.5f %4d %4d %4d\n";
    for $dim ("X", "Y", "Z") {
	$box->{$dim}{"LEN"} = $box->{$dim}{"hi"} - $box->{$dim}{"lo"};
    }

    print $datFile "Atoms\n\n";
    $index = 1;
    for $counter (sort numerically keys %{ $atoms } ) {
	%IMAGE = ();
	for $dim ("X", "Y", "Z") {
	    #if ($atoms->{$counter}{$dim . "COORD"} > $box->{$dim}{"hi"}) {
		#$IMAGE{$dim} = int($atoms->{$counter}{$dim . "COORD"}/$box->{$dim}{"LEN"}) + 1;
		#$IMAGE{$dim} = getImage($atoms->{$counter}{$dim . "COORD"}, $box->{$dim}{"hi"}, $box->{$dim}{"LEN"}, 1);
		#$atoms->{$counter}{$dim . "COORD"} -= (($IMAGE{$dim} - 1) * $box->{$dim}{"LEN"});
		#$atoms->{$counter}{$dim . "COORD"} -= ($IMAGE{$dim} * $box->{$dim}{"LEN"});
	    #} elsif ($atoms->{$counter}{$dim . "COORD"} < $box->{$dim}{"lo"}) {
		#$IMAGE{$dim} = -1 * int(abs($atoms->{$counter}{$dim . "COORD"})/$box->{$dim}{"LEN"}) - 1;
		#$IMAGE{$dim} = getImage($atoms->{$counter}{$dim . "COORD"}, $box->{$dim}{"lo"}, $box->{$dim}{"LEN"}, 0);
		#$atoms->{$counter}{$dim . "COORD"} -= ($IMAGE{$dim} * $box->{$dim}{"LEN"});
	    #} else {
		$IMAGE{$dim} = 0;
	    #}
	}
	#$IMAGE{X} = $IMAGE{Y} = $IMAGE{Z} = 0 if ($PARMS->{PARMS}{isTriclinic});

	$atoms->{$counter}{MOLECULEID} = $atoms->{$counter}{RESNUM} if ($ffType == 3); # mesodna fix
	$atoms->{$counter}{MOLECULEID} = $atoms->{$counter}{RESNUM};
	$out_string = sprintf($fmt, $index, $atoms->{$counter}{MOLECULEID}, $atoms->{$counter}{"PARMS"}{"INDEX"}, 
			      $atoms->{$counter}{"CHARGE"}, $atoms->{$counter}{"XCOORD"}, $atoms->{$counter}{"YCOORD"}, 
			      $atoms->{$counter}{"ZCOORD"}, $IMAGE{"X"}, $IMAGE{"Y"}, $IMAGE{"Z"});
	print $datFile $out_string;
	$index++;
    }

    print $datFile "\n";

}

sub getValenceParms {
    my ($atoms, $cons, $parms, $atmList) = @_;
    my ($i, $j, $k, $l, @currIndices, $PLIST, $tmp);

    for $i (@{ $atmList }) {
	if ($#{ $cons->{$i} } == 2) { #Inversion
	    @currIndices = ($i,@{ $cons->{$i} });
	    searchForValence($atoms, $parms->{INVERSIONS}, \@currIndices, "INVERSIONS", 1);
	}
	for $j (@{ $cons->{$i} }) {
	    next if ($j == $i);
	    if ($j > $i) {
		@currIndices = ($i, $j);
		searchForValence($atoms, $parms->{BONDS}, \@currIndices, "BONDS", 0);
	    }
	    for $k (@{ $cons->{$j} }) {
		next if ($k == $i || $k == $j);
		if ("${i}${j}${k}" > "${k}${j}${i}") {
		    @currIndices = ($i, $j, $k);
		    searchForValence($atoms, $parms->{ANGLES}, \@currIndices, "ANGLES", 0);
		}
		for $l (@{ $cons->{$k} }) {
		    next if ($l == $i || $l == $j || $l == $k);
		    if ("${i}${j}${k}${l}" > "${l}${k}${j}${i}") {
			@currIndices = ($i, $j, $k, $l);
			searchForValence($atoms, $parms->{TORSIONS}, \@currIndices, "TORSIONS", 0);
		    }
		}
	    }
	}
    }
}

sub searchForValence {
    my ($atoms, $parmList, $indices, $TYPE, $isInversion) = @_;
    my (@currTypes, $result, $error_code, $curr, $tmp); 
    my ($count, $IndexList, $valType, $bestParm, $i);

    @currTypes = ();
    $result = ();
    $error_code = "";
    for $count (@{ $indices }) {
	push @currTypes, $atoms->{$count}{FFTYPE};
	$error_code .= $atoms->{$count}{FFTYPE} . "-";
    }
    chop $error_code;
    if (! $isInversion) {
	@{ $IndexList->[0] } = @currTypes;
	@currTypes = reverse @currTypes;
	@{ $IndexList->[1] } = @currTypes;
    } else {
	$IndexList = [getPermutations(\@currTypes)];
    }

    for $curr (@{ $IndexList }) {
	($tmp, $count) = findValenceType($parmList, $curr, 0);
	if (keys %{ $tmp }) {
	    for $i (keys %{ $tmp }) {
		push @{ $result }, $tmp->{$i};
	    }
	}
    }

    if ($#{ $result } > -1) {
	$valType = '$' . $TYPE;
	$bestParm = getLeastX($result);
	for $i (values %{ $bestParm }) {
	    if (! $isInversion) { 
		next if exists($i->{ISCHILD});
		saveValence(\%{ eval($valType) }, $i, $indices);
	    } else { #have to figure out atom sequence for inversion only
		$IndexList = getCorrectInversionIndices($i, $indices, $atoms);
		for (@{ $IndexList }) {
		    saveValence(\%{ eval($valType) }, $i, $_);
		}
	    }
	}
	if ($TYPE eq "TORSIONS") {
	    ($indices->[1], $indices->[2]) = ($indices->[2], $indices->[1]) 
		if ($indices->[1] > $indices->[2]);
	    $ScaleTorsionList->{$indices->[1]}{$indices->[2]} = 0 
		if (! exists($ScaleTorsionList->{$indices->[1]}{$indices->[2]}));
	    $ScaleTorsionList->{$indices->[1]}{$indices->[2]}++;
	}
    }
    $ERRORS{$TYPE}{$error_code}++ if ($#{ $result } == -1);	
}

sub findValenceType {
    my ($MPARM, $type_keys, $count) = @_;
    my ($curr_key, $results, $curr_parm, %DAT, @tmp);
    my ($i, $xCount, $minXCount, @junk, $j, @new_keys);

    $curr_parm = $MPARM;
    $curr_key = $type_keys->[0];
    if (exists($curr_parm->{$curr_key}) || exists($curr_parm->{X})) {
        for $j ($curr_key,"X") {
            $curr_parm = $MPARM;
            next if (! exists($curr_parm->{$j}));
            $curr_parm = \%{ $curr_parm->{$j} };
            if ($#{ $type_keys } == 0) {
		$count++;
		$DAT{$count} = $curr_parm;
		return (\%DAT, $count);
            }
            @new_keys = @{ $type_keys };
            shift @new_keys;
            ($results, $count) = findValenceType($curr_parm, \@new_keys, $count);
            @tmp = keys %{ $results };
            if ($#tmp > -1) { # multiple so get one with least amount of Xs
		for $i (@tmp) {
		    $DAT{$i} = $results->{$i};
		}
            } 
        }
        return (\%DAT, $count);
    } else {
        return (undef, $count);
    }
}

sub getLeastX {
    my ($parmList) = $_[0];
    my ($i, @tmp, $xCount, $minXCount, $retParm, @pKeys, $j, $index);

    $minXCount = 99999;
    for $index (@{ $parmList }) {
	@tmp = keys %{ $index };
	next if (! @tmp);
	$i = $index->{$tmp[0]};
	@pKeys = split /\s+/, substr($i->{KEY},0,-1);
	$xCount = 0;
	for $j (@pKeys) {
	    $xCount++ if ($j eq "X");
	}
	if ($xCount < $minXCount) {
	    $retParm = $index;
	    $minXCount = $xCount;
	}
    }

    return $retParm;
}
sub saveValence {
    my ($TYPE, $VAL, $atomList) = @_;
    my ($rec, $i, $type, $counter, $j);
                                                                                                                                        
    $counter = $TYPE->{counter};
    $VAL->{USED} = 1;
    $rec->{DATA} = $VAL;
    for $i (@{ $atomList }) {
        push @{ $rec->{ATOMS} }, $i;
        $type = $ATOMS->{$i}{FFTYPE};
    }
    $counter++;
    $TYPE->{LIST}{ $counter } = $rec;
    if (exists($VAL->{NEXT}) and defined($VAL->{NEXT})) {
        for $i (@{ $VAL->{NEXT} }) {
            $rec = ();
            $i->{USED} = 1;
            $rec->{DATA} = $i;
            for $j (@{ $atomList }) {
                push @{ $rec->{ATOMS} }, $j;
            }
            $counter++;
            $TYPE->{LIST}{ $counter } = $rec;
        }
    }
    $TYPE->{counter} = $counter;
}
sub getCorrectInversionIndices {
    my ($inversionData, $atomList, $atomData) = @_;
    my ($centralIndex, $centralAtom, @inversionList, @indexCombo, $i, $j, @inversionKey, $sortedList);

    @inversionKey = split /\s+/, substr($inversionData->{KEY},0,-1);
    # for inversion of IJKL
    if ($inversionData->{TYPE} eq "UMBRELLA") { # umbrella torsion
	$centralIndex = 0; #central atom is index[1] J
    } elsif ($inversionData->{TYPE} eq "IT_JIKL") { # amber torsion
	$centralIndex = 2; #central atom is index[1] J
    } else { # charmm torsion
	$centralIndex = 0; # central atom is index[0] I
    }
    $centralAtom = $atomList->[0];
    shift @{ $atomList }; # remove the central atom index from the atom list
    shift @inversionKey; # remove the central atom type form the key list
    @{ $sortedList } = sort numerically @{ $atomList }; # sort the remainding atom indices in decreasing order
    $atomList = ();

    for $i (0 .. $#inversionKey) { # search of first atom of current atom type
	next if ($inversionKey[$i] eq "X");
	for $j (0 .. $#{ $sortedList }) { # atom with lower id checked first
	    if ($atomData->{ $sortedList->[$j] }{FFTYPE} eq $inversionKey[$i]) { # found atom type for current type
		$atomList->[$i] = $sortedList->[$j]; # save atom index in list
		splice @{ $sortedList }, $j, 1; # remove atom index from sorted list
		last; # exit sorted list loop
	    }
	}
    }

    for $i (0 .. $#inversionKey) {
	if (! $atomList->[$i]) { # atom key was not found, must have been an X
	    $atomList->[$i] = shift @{ $sortedList }; # place the smallest atom index
	}
    }

#    @{ $atomList } = sort numerically @{ $atomList }; 
    @indexCombo = ([$atomList->[0], $atomList->[1], $atomList->[2]],
		   [$atomList->[2], $atomList->[0], $atomList->[1]],
		   [$atomList->[1], $atomList->[2], $atomList->[0]]); # list for all 3 inversion about central atom
    for $j (0 .. $#indexCombo) {
	$atomList = $indexCombo[$j];
	for $i (0 .. ($centralIndex - 1)) { # add all atoms before central
	    $inversionList[$j][$i] = shift(@{ $atomList });
	}
	$inversionList[$j][$centralIndex] = $centralAtom; # add central atom
	while (@{ $atomList }) { # add all atoms after central
	    push @{ $inversionList[$j] }, shift @{ $atomList };
	}
	last if ($PARMS->{PARMS}{single_inversion}); # end if only using single inversion
    }

    return \@inversionList;
}

sub sortParmByIndex {
    my ($parm, $jVal) = @_;
    my (%INDICES, $i, $PLIST, $count, $j, $tmp);

    $PLIST = ();
    $count = 0;
    $tmp = 0;
    getValParent($parm, \%{ $PLIST }, \$tmp);
    $count = 0;
    for $i (keys %{ $PLIST }) {
	for $j (keys %{ $PLIST->{$i} }) {
	    next if (defined($jVal) and ($j != $jVal));
	    if (exists($PLIST->{$i}{$j}{INDEX})) {
		$count++;
		$INDICES{$PLIST->{$i}{$j}{INDEX}} = \%{ $PLIST->{$i}{$j} };
	    }
	}
    }
    return (\%INDICES, $count);
}

sub updateParmIndex {
    my ($parms) = $_[0];
    my ($i, $count, $PLIST, $index, @tmp, $j, $l, $k, $lammpsType);

    $parms->{PARMS}{USE_HBOND} = 0;
    $index = 1;
    for $i (keys %{ $parms->{ATOMTYPES} }) {
	$j = $parms->{ATOMTYPES}{$i}{TYPEID};
	$tmp[$j - 1] = $i;
    }
    $i = 0;
    while ($i <= $#tmp) {
	if (! $tmp[$i]) {
	    splice @tmp, $i, 1;
	    next;
	}
	$i++;
    }
    for $j (@tmp) {
	if (! $parms->{ATOMTYPES}{$j}{USED}) {
	    delete $parms->{ATOMTYPES}{$j};
	    delete $parms->{VDW}{$j};
	    for $i (@tmp) {
		next if (! exists($parms->{VDW}{$i}));
		delete $parms->{VDW}{$i}{$j} if (exists($parms->{VDW}{$i}{$j}));
	    }
	} else {
            $parms->{ATOMTYPES}{$j}{INDEX} = $index;
            $parms->{ATOMTYPES}{$j}{NAME} = $j;
	    $index++;
	}
    }

    for $j (@tmp) {
	    next if (! exists($parms->{VDW}{$j}));	
	    for $i (keys %{ $parms->{VDW}{$j} }) {
		next if (! exists($parms->{VDW}{$j}{$i}));
		for $l (keys %{ $parms->{VDW}{$j}{$i} }) {
                    $k = $parms->{VDW}{$j}{$i}{$l};
		    next if (! keys %{ $k });
                    if ($k->{TYPE} eq "DREIDHB") { # dreiding hb fix
                        unshift @{ $k->{VALS} }, $j;
                    }
		    if ($k->{TYPE} eq "LJ_6_12" and exists($parms->{PARMS}{is_tip4p})) { # tip4p fix
			$k->{Lammps}{name} = "lj/cut/coul/long/tip4p";
			$k->{Lammps}{opts} = "";
		    } 
		    if ($k->{IT} eq "hbond") {
			if (! $parms->{ATOMTYPES}{ $k->{VALS}[0] }{USED}) {
			    delete $parms->{VDW}{$j}{$i}{$l};
			    next;
			}
			if (! $parms->{BONDS}{$i}{$k->{VALS}[0]}{1}{USED} &&
			    ! $parms->{BONDS}{$k->{VALS}[0]}{$i}{1}{USED} &&
			    ! $parms->{BONDS}{$j}{$k->{VALS}[0]}{1}{USED} &&
                            ! $parms->{BONDS}{$k->{VALS}[0]}{$j}{1}{USED}) {
			    delete $parms->{VDW}{$j}{$i}{$l};
			    next;
			}
                        $k->{VALS}[0] = $parms->{ATOMTYPES}{ $k->{VALS}[0] }{INDEX}; # replace the hydrogen atom type with the type index
			if ($parms->{ATOMTYPES}{$j}{INDEX} <= $parms->{ATOMTYPES}{$i}{INDEX}) {
			    unshift @{ $k->{VALS} }, "i"; # specify that the donor is the type I
			} else {
			    unshift @{ $k->{VALS} }, "j"; # else the donor is type J
			}
			# format is i_type j_type hydrogen_type is_itype_donor d0 (alpha for morse) r0 cos_power 
			($k->{VALS}[0], $k->{VALS}[1]) = ($k->{VALS}[1], $k->{VALS}[0]);
            		$parms->{PARMS}{HAS_HBONDS} = 1;
            		$parms->{PARMS}{USE_HBOND} = 1;
			$parms->{PARMS}{USE_HBOND} = 2 if ($k->{TYPE} eq "MORSE_COSP");
                    }
		    $lammpsType = $k->{Lammps}{name};
		    $parms->{VDW}{TYPE}{$lammpsType} = $k->{Lammps}{opts};
                }
		next if ($i ne $j);
		for $k (values %{ $parms->{VDW}{$j}{$i} }) {
		    $k->{INDEX} = $parms->{ATOMTYPES}{$j}{INDEX};
		    $k->{DATA} = \%{ $parms->{ATOMTYPES}{$j} };
		}
	    }
    }
    $parms->{ATOMTYPES}{counter} = $parms->{VDW}{counter} = $index - 1;

    for $i ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
	$index = 0;
	$count = 0;
	$PLIST = ();
	getValParent($parms->{$i}, \%{ $PLIST }, \$count);
	for $j (sort numerically keys %{ $PLIST }) {
	    for $l (sort {$a cmp $b } keys %{ $PLIST->{$j} }) {
		if (! $PLIST->{$j}{$l}{USED} ) {
		    delete $PLIST->{$j}{$l};
		} else {
		    $index++;
		    $PLIST->{$j}{$l}{INDEX} = $index;
		    $lammpsType = $PLIST->{$j}{$l}{Lammps}{name} . " " . $PLIST->{$j}{$l}{Lammps}{opts};
		    $parms->{$i}{TYPE}{$lammpsType} = 1;
		}
	    }
	}
	removeEmptyParms($parms->{$i});
	$parms->{$i}{counter} = $index;
    }
}

sub removeEmptyParms {
    my ($parmData) = $_[0];
    my ($i);

    for $i (keys %{ $parmData }) {
	last if (exists($parmData->{$i}{INDEX}));
	next if ($i eq "counter" || $i eq "TYPE" || $i eq "Lammps");
	if (! keys %{ $parmData->{$i} }) {
	    delete $parmData->{$i};
	} else {
	    removeEmptyParms($parmData->{$i});
	    delete $parmData->{$i} if (! keys %{ $parmData->{$i} });
	}
    }
}

sub getPermutations {
    my ($inArray) = @_;
    my (@PERMS, $firstAtm, $i);

    $firstAtm = shift @{ $inArray };
    @PERMS = Permutate($inArray, []);
    for $i (0 .. $#PERMS) {
	unshift @{ $PERMS[$i] }, $firstAtm;
    }
    return @PERMS;
}

sub updateTorsionList {
    my ($torsionList) = $_[0];
    my ($i, $j, $k, @tmp, $l);
    my ($TLIST, $count);


    $count = 0;
    $TLIST = ();
    &getValParent($torsionList, \%{ $TLIST }, \$count);

    for $i (keys %{ $TLIST }) {
	for $j (keys %{ $TLIST->{$i} }) {
	    if ($j eq "") {
		delete $TLIST->{$i}{$j};
		$count--;
		next;
	    }	    
	    if ($TLIST->{$i}{$j}{NUM} > 1) { #multiple torsions
		@tmp = @{ $TLIST->{$i}{$j}{VALS} };
		$TLIST->{$i}{$j}{VALS} = ();
		for $l (0 .. ($TLIST->{$i}{$j}{PER} - 1)) {
		    $TLIST->{$i}{$j}{VALS}[$l] = $tmp[$l];
		}
		if ($TLIST->{$i}{$j}{Lammps}{name} eq "charmm") {
		    push @{ $TLIST->{$i}{$j}{VALS} }, getTorsionWeightFactor($TLIST->{$i}{$j});
		}
		$l = $TLIST->{$i}{$j}{PER};
		$k = 1;
		while ($k < $TLIST->{$i}{$j}{NUM}) {
		    %{ $TLIST->{$i}{"${j}${k}"} } = %{ $TLIST->{$i}{$j} };
		    $TLIST->{$i}{"${j}${k}"}{ISCHILD} = 1;
		    $TLIST->{$i}{"${j}${k}"}{INDEX} += $k;
		    $TLIST->{$i}{"${j}${k}"}{KEY} = "$k" . $TLIST->{$i}{"${j}${k}"}{KEY};
		    $TLIST->{$i}{"${j}${k}"}{VALS} = ();
		    for (0 .. ($TLIST->{$i}{$j}{PER} - 1)) {
			$TLIST->{$i}{"${j}${k}"}{VALS}[$_] = $tmp[$_ + $l];
		    }
		    if ($TLIST->{$i}{$j}{Lammps}{name} eq "charmm") {
			push @{ $TLIST->{$i}{"${j}${k}"}{VALS} }, 0;
		    }
		    $l += $TLIST->{$i}{$j}{PER};
		    push @{ $TLIST->{$i}{$j}{NEXT} }, \%{ $TLIST->{$i}{"${j}${k}"} };
		    $k++;
		}
	    } else {
		if ($TLIST->{$i}{$j}{Lammps}{name} eq "charmm") {
		    push @{ $TLIST->{$i}{$j}{VALS} }, getTorsionWeightFactor($TLIST->{$i}{$j});
		}
	    }
	}
    }
}

sub getTorsionWeightFactor {
    my ($torsion) = $_[0];
    my (@tmp, $atom1, $i, $j, $ts);

    if ($#{ $FILES } > 0 and $ffType > 4 and $ffType != 6) {
	return $torsion->{"1_4scale"};
    } elsif ($ffType == 1) { #amber
	return 0;
    } elsif ($ffType == 6) { #opls
	return 0;
    } elsif ($ffType == 2) { #charmm
	@tmp = split /\s+/, substr($torsion->{KEY},0,-1);
	for $i (0 .. $#tmp) {
	    next if ($tmp[$i] eq "X");
	    $atom1 = $tmp[$i];
	    $j = $i;
	    last;
	}
	if (! defined($atom1)) {
	    return 1;
	}
	return is5MemberRing($atom1, $j, \@tmp);
    } elsif ($ffType == 4) { #dreiding
	#return $PARMS->{PARMS}{"EXO_CYCLIC"}*$torsion->{"1_4scale"} 
	    #if (is5MemberRing($atom1, $j, \@tmp));
    } else {
	return 0;
    }
}

sub getTorsionScalingFactor {
    my ($torsions) = $_[0];
    my ($i, $j, $count, $atom1, $atom4, $index);
    $index = 0;
    for $i (1 .. $torsions->{counter}) {
	next if (! $torsions->{LIST}{$i}{DATA}{do_scale} or $torsions->{LIST}{$i}{DATA}{NUM} > 1);
	if (! $index) {
	    print "scaling torsions...";
	    $index = 1;
	}
	$atom1 = $torsions->{LIST}{$i}{ATOMS}[1];
	$atom4 = $torsions->{LIST}{$i}{ATOMS}[2];
	($atom1, $atom4) = ($atom4, $atom1) if ($atom1 > $atom4);
	$count = 1;
	$count = $ScaleTorsionList->{$atom1}{$atom4} if (exists($ScaleTorsionList->{$atom1}{$atom4}));
	#for $j (1 .. $torsions->{counter}) {
	    #$count++ if (($torsions->{LIST}{$j}{ATOMS}[1] == $atom1 && 
			  #$torsions->{LIST}{$j}{ATOMS}[2] == $atom4) ||
			 #($torsions->{LIST}{$j}{ATOMS}[2] == $atom1 &&
			  #$torsions->{LIST}{$j}{ATOMS}[1] == $atom4));
	#}

	next if ($count == 1);
	if (! exists($torsions->{LIST}{$i}{DATA}{scaled})) {
	    $torsions->{LIST}{$i}{DATA}{VALS}[0] /= $count;
	    $torsions->{LIST}{$i}{DATA}{scaled} = $count;
	} elsif ($torsions->{LIST}{$i}{DATA}{scaled} != $count) { #create new torsion
	    $torsions->{LIST}{$i}{DATA} = getNewParm($torsions->{LIST}{$i}{DATA}, "TORSIONS", $count, $i);
	}
    }
}

sub getNewParm {
    my ($parmData, $parmType, $count, $newKey) = @_;
    my ($newParm, $i, $DATA, $j, $newVal, @pkeys, $index);
    
    ($DATA, $j) = sortParmByIndex($PARMS->{$parmType});
    for $i (1 .. $PARMS->{$parmType}{counter}) {
	if ($DATA->{$i}{KEY} eq $parmData->{KEY} &&
	    exists($DATA->{$i}{scaled}) && 
	    $DATA->{$i}{scaled} == $count) {
	    return $DATA->{$i}; # found similar parm so return it
	}
    }
    
    #else create a new parm
    $newVal = $parmData->{VALS}[0] * $parmData->{scaled}/$count;
    @pkeys = split /\s+/, substr($parmData->{KEY},0,-1); 
    $index = $PARMS->{$parmType}{counter} + 1;
    $newParm = $PARMS->{$parmType}{shift @pkeys};
    for $i (0 .. ($#pkeys -1)) {
	$newParm = \%{ $newParm->{$pkeys[$i]} };
    }
    $newKey = pop(@pkeys) . $newKey;
    $newParm = \%{ $newParm->{$newKey} };
    $newParm = \%{ $newParm->{1} };
    %{ $newParm } = %{ $parmData };
    $newParm->{INDEX} = $index;
    $newParm->{do_scale} = 1;
    $newParm->{VALS} = ();
    @{ $newParm->{VALS} } = @{ $parmData->{VALS} };
    $newParm->{VALS}[0] = $newVal;
    $newParm->{scaled} = $count;
    $PARMS->{$parmType}{counter} = $index;

    return $newParm; #return the new parm
}

sub is5MemberRing {
    my ($atomType, $atomPos, $types) = @_;
    my ($atom1List, $l, $atom4List, $i, $j, $is5Member, $k);

    if ($atomPos == 0) { #already the first atom
	$atom1List = getAtmList($atomType, 0);
	$atom4List = getAtmList($atomType, 4);
    } elsif ($atomPos == 3) { #already the last atom
	$atom4List = getAtmList($atomType, 0);
	$atom1List = getAtmList($atomType, 4);
    } else {
	$atom1List = getAtmList($atomType, $atomPos);
	$atom4List = getAtmList($atomType, (3 - $atomPos));
    }

    $is5Member = 1;
  MAIN: for $i (values %{ $atom1List }) {
      next if ($ATOMS->{$i}{FFTYPE} ne $types->[0] && $types->[0] ne "X");
      for $j (values %{ $atom4List }) {
	  next if ($j == $i || ($ATOMS->{$j}{FFTYPE} ne $types->[3] && $types->[3] ne "X"));
	  for $k (@{ $BONDS->{$i} }) {
	      for $l (@{ $BONDS->{$k} }) {
		  if ($l == $i) {
		      $is5Member = 0;
		      last MAIN;
		  }
	      }
	  }
	  last MAIN;
      }
  }

    return $is5Member;
}

sub getAtmList {
    my ($atomType, $bondNum) = @_;
    my ($atom, $i, $BONDLIST, $count);

    for $i (keys %{ $ATOMS }) {
	next if ($ATOMS->{$i}{FFTYPE} ne $atomType);
	$atom = $i;
	last;
    }
    return undef if (! defined($atom));
    $BONDLIST = ();
    $count = 0;
    &getBondList(\%{ $BONDLIST }, $atom, \$bondNum, \$count);
    return $BONDLIST;
}

sub getBondList {
    my ($VALS, $atom, $bondC, $count) = @_;
    my ($i);

    if (${ $bondC } == 0) {
	${ $count }++;
	$VALS->{${ $count }} = $atom;
    } else {
	${ $bondC }--;
	for $i (@{ $CONS->{$atom} }) {
	    &getBondList($VALS, $i, $bondC, $count);
	}
    }
}

sub getValParent {
    my ($valList, $VList, $counter) = @_;
    my ($rec, $i, $j, $valid);

    for $i (keys %{ $valList }) {
	next if ($i eq "counter" || $i eq "TYPE");
	if (keys %{ $valList->{$i} }) {
	    $valid = 0;
	    for $j (keys %{ $valList->{$i} }) {
		next if ($i eq "counter" || $i eq "TYPE");
		if (exists($valList->{$i}{$j}{VALS})) {
		    $valid = 1;
		    last;
		}
	    }
	    if (! $valid) {
		getValParent($valList->{$i}, $VList, $counter);
	    } else {
		${ $counter }++;
		$VList->{${ $counter }} = \%{ $valList->{$i} };
	    }
	}
    }
}

sub printErrors {
    my ($errorlist, $fatal) = @_;
    my ($i, $j);

    for $i (keys %{ $errorlist }) {
	if ($fatal) {
	    print "\n---===$i ERRORS===----\n";
	} else {
	    print "\n---===MISSING $i TERMS===----\n";
	}
	for $j (keys %{ $errorlist->{$i} }) {
	    print "$j: (occurred $errorlist->{$i}{$j} times)\n";
	}
    }

    die "The script cannot contine\n" if ($fatal);
}

sub addLammpsParms {
    my ($parms, $atoms, $mols, $angles) = @_;
    my ($i, $parmHybrid, %offDiag, @tmp, $k, $j, $index1, $index2, $hybridOpt, $l);

    $parms->{PARMS}{SOLUTE} = GetSoluteAtoms($atoms, $mols);
    $parms->{PARMS}{SUFFIX} = $suffix;
    $parms->{PARMS}{FFTYPE} = $ffType;
    $parms->{PARMS}{NUM_FILES} = $#{ $FILES } + 1;
    @tmp = keys %{ $inputType->{names} };
    $parms->{PARMS}{INPUTNAME} = "@tmp";
    $parms->{PARMS}{INPUTLOC} = $inputType->{loc};
    $parms->{PARMS}{NUM_ATOMS} = scalar(keys %{ $atoms });

    for $i (1 ... $angles->{counter}) {
	if ($#{ $ANGLES->{LIST}{$i}{ATOMS} } == 2) {
	    if ($atoms->{ $angles->{LIST}{$i}{ATOMS}[0] }{ELENUM} == 1 &&
		$atoms->{ $angles->{LIST}{$i}{ATOMS}[1] }{ELENUM} == 8 &&
		$atoms->{ $angles->{LIST}{$i}{ATOMS}[2] }{ELENUM} == 1) {
		$parms->{PARMS}{SHAKE_ANGLE} = " a $angles->{LIST}{$i}{DATA}{INDEX}";
		$parms->{PARMS}{SHAKE_MASS} = $parms->{ATOMTYPES}{ $atoms->{ $angles->{LIST}{$i}{ATOMS}[0] }{FFTYPE} }{MASS};
		last;
	    }
	}
    }    
    if (! exists($parms->{PARMS}{SHAKE_MASS})) {
	for $i (keys %{ $atoms }) {
	    if ($atoms->{$i}{ELENUM} == 1) {
		$parms->{PARMS}{SHAKE_MASS} = $parms->{ATOMTYPES}{ $atoms->{$i}{FFTYPE} }{MASS};
	  	last;
	    }
	}
    }
    @tmp = keys %{ $parms->{ATOMTYPES} };
    $parmHybrid = determineIfHybrid($parms->{VDW});
    $parms->{PARMS}{HYBRID_OVERLAP} = 1 if ($parmHybrid eq "hybrid/overlap ");
    for $i (keys %{ $parms->{VDW} }) {
	next if ($i eq "TYPE" || $i eq "counter");
	for $k (keys %{ $parms->{VDW} }) {
	    #next if (($k eq $i && scalar keys %{ $parms->{VDW}{$i}{$k} } == 1) || $k eq "TYPE" || $k eq "counter");
	    next if ($k =~ /type|counter/i);
	    next if (! exists($parms->{VDW}{$i}{$k}));
	    next if ($k eq $i && ! $parmHybrid);
	    @tmp = sort numerically keys  %{ $parms->{VDW}{$i}{$k} };
	    for $l (@tmp) {
		if ($parms->{ATOMTYPES}{$i}{INDEX} < $parms->{ATOMTYPES}{$k}{INDEX}) {
		    $index1 = $parms->{ATOMTYPES}{$i}{INDEX};
		    $index2 = $parms->{ATOMTYPES}{$k}{INDEX};
		} else {
		    $index2 = $parms->{ATOMTYPES}{$i}{INDEX};
		    $index1 = $parms->{ATOMTYPES}{$k}{INDEX};
		}
		if ($parmHybrid) {
		    $hybridOpt = sprintf("%-18s",$parms->{VDW}{$i}{$k}{$l}{Lammps}{name});
		} else {
		    $hybridOpt = "";
		}
		$offDiag{$index1}{$index2} .= sprintf("%-15s %-4d %-4d $hybridOpt","pair_coeff",$index1,$index2);
		for $j (@{ $parms->{VDW}{$i}{$k}{$l}{VALS} }) {
		    if (IsDecimal($j)) {
			$offDiag{$index1}{$index2} .= sprintf("%25.15f ", $j);
		    } elsif($j =~ /[a-zA-Z]/) {
			$offDiag{$index1}{$index2} .= sprintf("%25s ", $j);
		    } else {
			$offDiag{$index1}{$index2} .= sprintf("%25d ", $j);
		    }
		}
		$offDiag{$index1}{$index2} .= "\n";
	    }
	}
    }
    for $i (sort numerically keys %offDiag) {
	for $j (sort numerically keys %{ $offDiag{$i} }) {
	    $parms->{PARMS}{OFF_DIAG} .= $offDiag{$i}{$j};
	}
    }

    #reax fix
    if ($ffType == 5) {
	$parms->{PARMS}{REAX} = 1;
	@tmp = ();
	for $i (keys %{ $parms->{ATOMTYPES} }) {
	    next if ($i eq "counter");
	    $j = $parms->{ATOMTYPES}{$i}{INDEX};
	    $tmp[$j - 1] = $parms->{ATOMTYPES}{$i}{TYPEID};
	}
	$j = "";
	for $i (@tmp) {
	    $j .= "$i ";
	}
	$parms->{PARMS}{OFF_DIAG} .= sprintf("%-15s %-4s %-4sffield.reax $j\n\n", "pair_coeff", "*", "*");
    }
}
	    	
sub checkAtomTypes {
    my ($atoms, $atomTypes, $ffTypeID) = @_;
    my ($i, $ELEMENTS, $ffType, $atmName, $eleList);
    
    $ELEMENTS = &LoadElements;
    for $i (keys %{ $ELEMENTS }) {
        $eleList->{uc($ELEMENTS->{$i}{SYMBOL})} = $i;
    }


    for $i (keys %{ $atoms }) {
	$totAtms++;
	$ffType = $atoms->{$i}{FFTYPE};
	$atmName = $atoms->{$i}{ATMNAME};
        if (exists($atomTypes->{$ffType})) {
            ($atoms->{$i}{ELENAME}, $atoms->{$i}{ELENUM}) = FindElement($ffType, $atomTypes, $ELEMENTS);
            $atomTypes->{$ffType}{USED} = 1;
            $atoms->{$i}{PARMS} = \%{ $atomTypes->{$ffType} };
            $atoms->{$i}{CHARGE} = $atomTypes->{$ffType}{CHARGE} if ($atomTypes->{$ffType}{USE_CHARGE});

	}elsif ($ffTypeID == 5) {
            ($ffType, $atoms->{$i}{ELENUM}) = getElement($eleList, $ffType, $atmName);
            if (! defined($ffType) or ! exists($atomTypes->{$ffType})) {
                $ERRORS{ELEMENT}{$atoms->{$i}{FFTYPE}}++;
                next;
            }
	    $atoms->{$i}{ELENAME} = $ffType;
            $atoms->{$i}{FFTYPE} = $ffType;
            $atomTypes->{$ffType}{USED} = 1;
            $atoms->{$i}{PARMS} = \%{ $atomTypes->{$ffType} };

	} else {
	    $ERRORS{ATOMTYPES}{$ffType}++;
	}
    }
}

sub getElement {
    my ($eleList, $ffType, $atmName) = @_;
    my ($i, $element);

    for $i ($ffType, $atmName) {
	$element = $i;
	next if ($element !~ /^([A-Za-z])([a-z]?)/);
	$element = $1;
	$element .= $2 if ($2);
	return ($element, $eleList->{uc($element)}) if (exists($eleList->{uc($element)}));
	$element = $1;
	return ($element, $eleList->{uc($element)}) if (exists($eleList->{uc($element)}));
    }

}

sub determineIfHybrid {
    my ($parms) = $_[0];
    my (@tmp, $i, $k, $j, $l, $jType, $kType, $count, %vdwType, @tmp1);
    my ($parmHybrid) = "";

    @tmp = grep {!/TYPE|counter/i} keys %{ $parms };
    for $i (@tmp) {
        for $k (keys %{ $parms->{$i} }) {
	    for $l (keys %{ $parms->{$i}{$k} }) {
		$vdwType{ $parms->{$i}{$k}{$l}{TYPE} } = 1;
	    }
	}
    }

    $parmHybrid = "hybrid " if (scalar keys %vdwType > 1);
    for $i (keys %vdwType) {
	last if (! $parmHybrid);
        for $j (@tmp) {
            $jType = findVDWEntry($parms->{$j}{$j}, $i);
            next if (! $jType);
            for $k (@tmp) {
		next if (exists($parms->{$j}{$k}));
                next if ($k eq $j or $PARMS->{ATOMTYPES}{$k}{INDEX} > $PARMS->{ATOMTYPES}{$j}{INDEX});
                $kType = findVDWEntry($parms->{$k}{$k}, $i);
                next if (! $kType or findVDWEntry($parms->{$j}{$k}, $i) or findVDWEntry($parms->{$k}{$j}, $i));
		$count = scalar(keys %{ $parms->{$j}{$k} }) + 1;
		$parms->{$j}{$k}{$count}{ATOM} = $parms->{$j}{$k}{$count}{KEY} = "${j} ${k} ";
		%{ $parms->{$j}{$k}{$count}{Lammps} } = %{ $parms->{$j}{$j}{$jType}{Lammps} };
		$parms->{$j}{$k}{$count}{USED} = 1;
		$parms->{$j}{$k}{$count}{TYPE} = $i;
		$parms->{$j}{$k}{$count}{VALS} = getPairMix($parms->{$j}{$j}{$jType}{VALS}, 
							    $parms->{$k}{$k}{$kType}{VALS}, $i);
            }
        }
    }

    for $i (@tmp) {
	last if ($PARMS->{PARMS}{HYBRID_OVERLAP});
	for $k (keys %{ $parms->{$i} }) {
	    @tmp1 = keys %{ $parms->{$i}{$k} };
	    if(scalar(@tmp1) > 1) {
		$parmHybrid = "hybrid/overlap ";
		last;
	    }
	    for $j (@tmp1) {
		if ($j > 1) {
		    $parmHybrid = "hybrid/overlap ";
		    last;
		}
            }
        }
    }

    return $parmHybrid;
}

sub getPairMix {
    my ($iVals, $jVals, $pairType) = @_;
    my (@VALS, $mixType);

    $mixType = $PARMS->{PARMS}{mix_rule};

    if ($mixType eq "geometric") {
	$VALS[0] = sqrt($iVals->[0]*$jVals->[0]);
    } else {
	$VALS[0] = 0.5*($iVals->[0]+$jVals->[0]);

    }
    $VALS[1] = sqrt($iVals->[1]*$jVals->[1]);
  
    if ($pairType ne "LJ_12_10" and $#{ $iVals } == 2) {
        $VALS[2] = 0.5*($iVals->[2]+$jVals->[2]);
    }
	
    return \@VALS;
}

sub findVDWEntry {
    my ($vdw, $lammpsName) = @_;
    my ($i, $returnVal);

    for $i (keys %{ $vdw }) {
	if ($vdw->{$i}{TYPE} eq $lammpsName) {
	    $returnVal = $i;
	    last;
	}
    }

    return $returnVal;
}

sub getImage {
    my ($atomPos, $boxDim, $boxLen, $isHi) = @_;
    my ($imageVal) = 0;

    if ($isHi) {
	while ($boxDim < $atomPos) {
	    $imageVal++;
	    $atomPos -= $boxLen;
	}
    } else {
	while ($atomPos < $boxDim) {
	    $imageVal--;
	    $atomPos += $boxLen;
	}
    }
    return $imageVal;
}

sub findDuplicateParms {
    my ($parms) = $_[0];
    my ($PLIST, $i, $searchStr, $count, $j, $k, $SLIST);

    for $i ("VDW", "BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
	$SLIST = ();
	$count = 0;
	$PLIST = ();
        getValParent($parms->{$i}, \%{ $PLIST }, \$count);
	for $j (keys %{ $PLIST }) {
	    for $k (keys %{ $PLIST->{$j} }) {
		$searchStr = "$PLIST->{$j}{$k}{KEY} $PLIST->{$j}{$k}{TYPE} @{$PLIST->{$j}{$k}{VALS}}";
		if (exists($SLIST->{$searchStr})) {
		    delete $PLIST->{$j}{$k};
		} else {
		    $SLIST->{$searchStr} = 1;
		}
	    }
	}
     }
}

sub createLammpsClusterScript {
    my ($ffs) = $_[0];
    my ($ffList, $scriptFile, $scriptCmd, $reaxff);
   
    if ($#{ $ffs } == 0) {
	$ffList = $ffs->[0]{FF};
	$reaxff = $ffs->[0]{FF} if ($ffs->[0]{FFTYPE} eq "REAX");
    } else {
	$ffList = '"';
	for (@{ $ffs }) {
	    $ffList .= "$_->{FF} ";
	}
	chomp $ffList;
	$ffList .= '"';
    }

    $scriptFile = "$Bin/createClusterScript.pl";
    if (! -e $scriptFile) {
	print "Failed!\n";
	return;
    }
    $scriptFile .= " -r $reaxFF" if (defined($reaxff));

    $scriptCmd = "$scriptFile -p $suffix -b $bgfFile -f $ffList -i in.${suffix} -d data.${suffix} -s ${suffix}.lammps.pbs >& _createscript";
    system($scriptCmd);
    system("rm -fr _createscript");
}

sub getTIP4Popts {
    my ($parms, $atoms, $bonds) = @_;
    my ($tip4Str, $oType, $hType, $bType, $aType);

    # first get the type of the ow/hw atoms
    for $i (keys %{ $parms->{ATOMTYPES} }) {
	if ($i =~ /OW/) {
	    $oType = $parms->{ATOMTYPES}{$i}{INDEX};
	} elsif ($i =~ /HW/i) {
	    $hType = $parms->{ATOMTYPES}{$i}{INDEX};
	}
     }

     #get ow - hw bond type
     for $i (1 ... $bonds->{counter}) {
         if ($bonds->{LIST}{$i}{DATA}{KEY} =~ /HW OW/ or $bonds->{LIST}{$i}{DATA}{KEY} =~ /OW HW/) {
  	     $bType = $bonds->{LIST}{$i}{DATA}{INDEX};
             last;
         }
     }
     die "ERROR: OW - HW bond type not found for TIP4P water!\n" if (! defined($bType));
     #angle type is already recorded for shake
     $aType = $parms->{PARMS}{SHAKE_ANGLE};
     $aType =~ /(\d+)/;
     $aType = $1;

     $tip4Str = "lj/cut/coul/long/tip4p $oType $hType $bType $aType $parms->{PARMS}{tip4_om_dist} 10";
     delete $parms->{VDW}{TYPE}{"lj/cut/coul/long/tip4p"};
     $parms->{VDW}{TYPE}{$tip4Str} = "";
}

sub getInputTypes {
    my ($isStr) = $_[0];

    my ($searchDir) = "$Bin/dat/LAMMPS/";
    my (@list) = `ls ${searchDir}/in.lammps.*`;
    my ($i, $typeStr);

    for $i (@list) {
        chop $i;
	if ($i =~ /in.lammps.(\S+)/) {
	    if ($isStr) {
		$typeStr .= "$1 ";
	    } else {
		$typeStr->{$1} = "$i";
	    }
	}
    }
    return $typeStr;
}

sub setOpts {
    my ($bgf, $parms, $opts) = @_;
    my ($i, $j, $hasCharge, $tmp, $count, $add);

    $parms->{PARMS}{OPTIONS}{PERIODICITY} = 3;
    $parms->{PARMS}{OPTIONS}{PERIODICITY} = 0 if ($opts =~ /finite/i);
    $parms->{PARMS}{OPTIONS}{PERIODICITY} = 2 if ($opts =~ /2d|slab|dimension 2/i);
    $parms->{PARMS}{OPTIONS}{SHAKE} = 1;
    $parms->{PARMS}{OPTIONS}{SHAKE} = 0 if ($opts =~ /no shake/i);
    $parms->{PARMS}{OPTIONS}{SHAKE_WHAT} = $1 if ($opts =~ /shake (solute|solvent)/i);
    $parms->{PARMS}{OPTIONS}{VDW_EWALD} = 0;
    $count = 0;
    getValParent($parms->{VDW}, \%{ $tmp }, \$count);

    $hasCharge = 0;
    for $i (values %{ $bgf }) {
        if ($i->{CHARGE} > 0 or $i->{CHARGE} < 0) {
            $hasCharge = 1;
            last;
        }
    }

    if ($opts =~ /vdw ewald|ewald vdw/i) {
	$add = "long long";
	$add = "long off" if (! $hasCharge);
	$parms->{PARMS}{OPTIONS}{VDW_EWALD} = 1;
        for $i (1 .. $count) {
	    for $j (values %{ $tmp->{$i} }) {
        	$j->{Lammps}{name} =~ s/lj\/charmm\/coul\/long\S*/lj\/coul $add/g;
        	$j->{Lammps}{name} =~ s/lj\/cut\/coul\/long\S*/lj\/coul $add/g;
        	$j->{Lammps}{name} =~ s/buck\/coul\S*/buck\/coul $add/g;
		$j->{Lammps}{opts} = " $parms->{PARMS}{cut_vdw} $parms->{PARMS}{cut_coul}";
	    }
	}
    }

    if (! $hasCharge) {
	for $i (1 .. $count) {
	    for $j (values %{ $tmp->{$i} }) {
		$j->{Lammps}{name} =~ s/\/coul\S+//;
		$j->{Lammps}{name} = "lj/gromacs" if ($j->{Lammps}{name} eq "lj/charmm");
	    }
 	}
    } elsif ($parms->{PARMS}{OPTIONS}{PERIODICITY} == 0) {
        for $i (1 .. $count) {
	    for $j (values %{ $tmp->{$i} }) {
		$j->{Lammps}{name} =~ s/lj\/charmm\/coul\/long\/opt/lj\/charmm\/coul\/charmm/g;
		$j->{Lammps}{name} =~ s/\/coul\/long\ /\/coul\/cut /g;
	    }
        }
    }	
}

sub init {
    my ($i, $FFILES, $inputFFType, $inputList, $inputStr, $count);
    
    getopt('bfstrmlo',\%OPTS);
    die &usage if (! exists($OPTS{b}) || ! exists($OPTS{f}));
    ($bgfFile, $FF, $suffix, $inputStr, $reaxFF, $inputFFType, $writeLabels, $opts) = 
	($OPTS{b},$OPTS{f},$OPTS{s},$OPTS{t},$OPTS{r},$OPTS{m},$OPTS{l},$OPTS{o});
    
    print "Step 1: Initializing...";
    FileTester($bgfFile);
    if (lc($FF) =~ /reax/ and ! defined($reaxFF)) {
	print "using default Reax force field...";
	$reaxFF = "/home/noische/ff/Reax.ff";
	#$inputStr = "reax";
    }
    #$FF = "$reaxFF" if (lc($FF) =~ /reax/); #for now only allow pure reax simulations
    $reaxFF =~ s/^\s+// if (defined($reaxFF));
    $reaxFF =~ s/\s+.*$// if (defined($reaxFF)); # only allow 1 reax ff to be specified

    ($FFILES, $ffType) = ReadFFs($FF, $reaxFF);
    $inputStr = "mesodna" if ($ffType == 3);
   
    $writeLabels = 1 if (! defined($writeLabels));
    if ($writeLabels !~ /0|no/) {
	$writeLabels = 1;
    } else {
	$writeLabels = 0;
    }
    $suffix = "lammps" if (! defined($suffix));
    $inputList = getInputTypes(0);
    $inputStr = "min heat md 2pt" if (! defined($inputStr)); 
    while($inputStr =~ /(\S+)/ig) {
	if(-e $1) {
	    $inputType->{names}{custom} = 1;
	    push @{ $inputType->{loc} }, $1
	} elsif (exists($inputList->{$1})) {
	    $inputType->{names}{$1} = 1;
	    push @{ $inputType->{loc} },$inputList->{$1};
	}
    }
    $BONDS->{counter} = $ANGLES->{counter} = $TORSIONS->{counter} = $INVERSIONS->{counter} = 0;
    $totAtms = 0;
    print "Done\n";
    $ffType = $inputFFType if (defined($inputFFType));
    return $FFILES;
}
    
sub usage {
    my ($list) = getInputTypes(1);
    
return <<DATA;
This script will generate LAMMPS data and input files from a bgf structure
usage: $0 -b bgfFile -f \"ff1 ff2...\" -s [suffix] -o [options] -t [sim type] -r [reaxff]
	-b bgfFile: location of BGF file. Must have correct atom types
	-f \"forcefield1 forcefield2...\": 1 or more Cerius2/Polygraf formatted forcefields
		Valid entries are
		AMBER91/96/99 - the AMBER 1991/1996/1999 forcefield for proteins and DNA
		AMBER03 - the AMBER 1999 (+ 2003 modifications) force field
		GAFF - the General AMBER forcefield for small molecules
		CHARMM - the CHARMM par_all27_prot_na Protein/Nucleic Acid forcefield
		CHARMM_LIPID - the CHARMM par_all27_prot_lipid Protein/Lipid force field
		DREIDING - The 1990 DREIDING forcefield with F3C waters
		MESODNA - The DNA Meso scale forcefield v.6.0
		REAXFF - The REAX force field. You will have to specify a location with -r
		--or-- you can specify your own forcefield location
		NOTE: you can specify multiple forcefields by enclosing them in ""
		NOTE: If specifying multiple forcefields with the same atom types,
			the atom type values will come from the last matching forcefield.
	-s [suffix]: (optional) When specified, the program will generate in.[suffix]
		and data.[suffix] as the files. If not specified, the output will be
		in.lammps, data.lammps
	-o [options]: (optional). Controls various options in input file. Valid entries include:
		"2D or dimension 2" - for 2D simulations
		"finite" - for 0D (isolated) simulations
		"no shake|shake solute/solvent" - shake constraints on the hydrogens are turned on by default. 
			This turns it off or applies it to the solute or solvent only
		"ewald vdw" - calculate long range vdw using ewald summations. Only valid for lj and exp6 potentials.
	-t [sim]: (optional). Specifies the type of input file to create. See $Bin/dat/LAMMPS for a list
		Current options include "$list"
		or you can specify your own input file

Report any bugs to tpascal\@wag.caltech.edu
DATA

}
