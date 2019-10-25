#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use Packages::General qw(FileTester LoadFFs LoadElements getFFtype);
use Packages::Math::Amoeba qw(MinimiseND);
use File::Basename qw(basename dirname);
use Getopt::Std qw(getopt);
use strict;
use Cwd 'abs_path';

sub updateAtomTypeData;
sub createLammpsFiles;
sub getOptimalParams;
sub calcResiduals;
sub getPairCoeffs;
sub mix;
sub getNewSearchDirection;
sub assignBGFS;
sub getObjectiveFunctionData;
sub init;
sub createFF;
sub offDiag;
sub getFileEngData;
sub getNewPairCoeffs;
sub getParmOrder;

my ($energyFile, $lammpsBinary, $saveName, $lammpsInput, $PARMS, $FREEZE, $optOffDiag, $DATFILE);
my ($prefix, $DATA, $ffData, $createLammpsScript, $optScale, $EPROFILE, $FILES, $count);

$|++;
&init;
$PARMS = LoadFFs($ffData, undef, 0);
$DATA = &getFileEngData($FILES);
print "Creating LAMMPS data files...";
&createLammpsFiles($DATA, $ffData);
&updateAtomTypeData($DATA, $PARMS);
&offDiag($PARMS) if ($optOffDiag);
print "Done\nCreating objective function...";
($PARMS->{OPT}{guess}, $PARMS->{OPT}{scale}, $PARMS->{OPT}{map}) = getBounds($PARMS->{VDW}, $FREEZE, $optScale);
$PARMS->{OBJ} = $DATA;
print "Done\n";
&getOptimalParams($PARMS, $saveName);
print "Creating optimized FF $saveName...";
&createFF($PARMS, $saveName, $ffData);
print "Done\n";

sub getFileEngData {
    my ($flist) = $_[0];
    my ($i, $bgfs, $bgfLoc, $findCmd, $engFile, $data, $index, $rec, %DATA, $j);

    $index = 0;
    for $i (@{ $flist }) {
	$bgfs = ();
	$data = ();
	$bgfLoc = $i->{BGF};
	$engFile = $i->{ENG};
	if (-d $bgfLoc) {
	    $findCmd = "find $bgfLoc -name '*.bgf' -print";
	} else {
	    my ($d, $n) = (dirname($bgfLoc), basename($bgfLoc));
	    $findCmd = "find $d -name '" . $n . "' -print";
	}
	if (! open(FINDCMD, "$findCmd |")) {
	    die "ERROR: No valid bgf files found while searching \"$bgfLoc\"!\n";
	}
	while (<FINDCMD>) {
	    chomp;
	    push @{ $bgfs }, $_;
	}
	close FINDCMD;
	print "Parsing Energy function file $engFile...";
	$data = getObjectiveFunctionData($engFile);
	&assignBGFS($data, $bgfs);
	for $j (keys %{ $data }) {
	    $rec = $data->{$j};
	    $DATA{$index} = $rec;
	    $index++;
	}
	print "Done\n";
    }

    return \%DATA;
}

sub offDiag {
    my ($parms) = $_[0];
    my ($i, $j, @vals, $aa, $bb);

    for $i (keys %{ $parms->{VDW} }) {
	for $j (keys %{ $parms->{VDW} }) {
	    next if ($i eq $j || $i lt $j);
	    next if (exists($parms->{VDW}{$i}{$j}) or exists($parms->{VDW}{$j}{$i}));
	    $aa = $parms->{VDW}{$i}{$i}{1};
	    $bb = $parms->{VDW}{$j}{$j}{1};
	    for (keys %{ $aa }) {
		next if ($_ eq "VALS");
		$parms->{VDW}{$i}{$j}{1}{$_} = $aa->{$_};
	    }
	    if (uc($aa->{TYPE}) eq "VDW_MORSE") {
		$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
		$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
		$vals[2] = 0.5*($aa->{VALS}[2]+$bb->{VALS}[2]);
	    } elsif (uc($aa->{TYPE}) eq "LJ_6_12") {
		$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
		$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
		if ($#{ $aa->{VALS} } > 1) {
		    $vals[2] = sqrt($aa->{VALS}[2]*$bb->{VALS}[2]);
		    $vals[3] = 0.5*($aa->{VALS}[3]+$bb->{VALS}[3]);
		}
	    } elsif (uc($aa->{TYPE}) =~ /EXPO_6|BUCKINGHAM/) {
		$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
		$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
		$vals[2] = sqrt($aa->{VALS}[2]*$bb->{VALS}[2]);
	    }

	    @{ $parms->{VDW}{$i}{$j}{1}{VALS} } = @vals;
	}
    }
}

sub getNewPairCoeffs {
    my ($data, $ff_type) = @_;
    my ($i, $j, $diag, $offDiag, $curr);

    for $i (keys %{ $data->{VDW} }) {
	for $j (keys %{ $data->{VDW}{$i} }) {
	    $curr = $data->{VDW}{$i}{$j}{1};
	    if ($ff_type eq "MPSIM") {
		$diag .= sprintf("%-4s   %3d", $i, $curr->{NUM}) if ($i eq $j);
		$offDiag .= sprintf("%-4s -%-4s   %3d", $i, $j, $curr->{NUM}) if ($i ne $j);
		for (@{ $curr->{VALS} }) {
		    $offDiag .= sprintf("%10.5f",$_) if ($i ne $j);
		    $diag .= sprintf("%10.5f",$_)  if ($i eq $j);
		}
		$diag .= "\n" if ($i eq $j);
		$offDiag .= "\n" if ($i ne $j);
	    } elsif ($ff_type eq "CERIUS2") {
		$diag .= sprintf(" %-12s%-13s",$i,$curr->{TYPE}) if ($i eq $j);
		$offDiag .= sprintf(" %-9s%-9s%-13s",$i,$j,$curr->{TYPE}) if ($i ne $j);
		for (@{ $curr->{VALS} }) {
		    $diag .= sprintf("%8.4G",$_) if ($i eq $j);
		    $offDiag .= sprintf("%8.4G",$_)  if ($i ne $j);
		}
		$diag .= "\n" if ($i eq $j);
		$offDiag .= "\n" if ($i ne $j);
	    }
	}
    }

    return ($diag, $offDiag);
}

sub createFF {
    my ($data, $file, $ffFiles) = @_;
    my ($ff_type, $offset, $outData, $buf, $valid, $curr, $i, $j, $tot, $delimiter, $tmp);
    my ($diag, $offDiag);

    $ff_type = getFFtype($ffData->[0]);
    $tot = (-s $ffData->[0]);
    if ($ff_type eq "CERIUS2") {
	$offset->{1}{START} = `grep -b DIAGONAL_VDW $ffFiles->[0]`;
	$offset->{2}{START} = `grep -b OFF_DIAGONAL_VDW $ffFiles->[0]`;
	$delimiter = "END";
    } elsif ($ff_type eq "MPSIM") {
	$offset->{1}{START} = `grep -b "VDW AT ITY" $ffFiles->[0]`;
	$offset->{2}{START} = `grep -b "NONBOND-OFF" $ffFiles->[0]`;
	$delimiter = "*";
    }
    $tmp = `grep -b '$delimiter' $ffFiles->[0]`;
    for $i (keys %{ $offset }) {
	$offset->{$i}{START} =~ /^(\d+)/;
	$offset->{$i}{START} = $1;
	while ($tmp =~ /(\d+)/g) {
	    if ($1 > $offset->{$i}{START}) {
		$offset->{$i}{STOP} = $1;
		last;
	    }
	}
    }

    open FF, $ffFiles->[0] or die "ERROR: Cannot open $ffFiles->[0]: $!\n";
    read FF, $buf, ($offset->{1}{START} - 1), 0;
    $outData = $buf;
    undef($buf);
    if ($ff_type eq "MPSIM") {
	$outData .= "\nVDW AT ITY       RNB      DENB     SCALE\n";
    } elsif ($ff_type eq "CERIUS2") {
	$outData .= "DIAGONAL_VDW\n";
    }
    ($diag, $offDiag) = getNewPairCoeffs($data, $ff_type, 1);
    $outData .= $diag;

    if ($offset->{2}{START}) {
	$buf = "";
	seek FF,$offset->{1}{STOP},0;
	read FF,$buf,($offset->{2}{START}-$offset->{1}{STOP}),0;
	$outData .= $buf;
	if ($ff_type eq "CERIUS2") {
	    $outData .= "\nOFF_DIAGONAL_VDW\n";
	} elsif ($ff_type eq "MPSIM") {
	    $outData .= "*\nNONBOND-OFF TYPE      RVDW      DVDW     SCALE\n";
	}
	$outData .= $offDiag if(defined($offDiag));
    }
    if ($ff_type eq "CERIU2") {
	$outData .= "END\n";
    } elsif ($ff_type eq "MPSIM") {
	$outData .= "*\n";
    }
    $buf = "";
    seek  FF,$offset->{2}{STOP},0;
    read FF,$buf,($tot-$offset->{2}{STOP}),0;
    $outData .= $buf;

    open OUTDATA, "> $file" or die "ERROR: Cannot write to $file: $!\n";
    print OUTDATA $outData;
    close OUTDATA;

    $file =~ s/\.\w+$//;
    system("mv _lmptmp/in.tmp ./in.${file}");
}

sub updateAtomTypeData {
    my ($data, $ffData) = @_;
    my ($i, @tmp, $datfile, $start, $j, $ffTypeList);
    
    $ffTypeList = ();

    @tmp = keys %{ $data };
    $datfile = $data->{ $tmp[0] }{LMPDAT};
    $start = 0;
    open DAT, $datfile or die "ERROR: Cannot open $datfile: $!\n";
    while (<DAT>) {
	chomp;
	$start = 1 if ($_ =~ /^Masses/);
	last if ($_ =~ /Coeff/);
	if ($start and $_ =~ /(\d+)\s+\d+\.\d+\s+\#\s+(\w+)/) {
	    $ffTypeList->{$2} = $1;
	}
    }
    close DAT;
    die "ERROR: Cannot parse LAMMPS data file $datfile!\n" if (! defined($ffTypeList));

    for $i (keys %{ $ffData->{ATOMTYPES} }) {
	if (! exists($ffTypeList->{$i})) {
	    delete $ffData->{ATOMTYPES}{$i};
	    delete $ffData->{VDW}{$i};
	    for $j (keys %{ $ffData->{VDW} }) {
		delete $ffData->{VDW}{$j}{$i} if exists($ffData->{VDW}{$j}{$i});
	    }
	} else {
	    $ffData->{ATOMTYPES}{$i}{INDEX} = $ffTypeList->{$i};
	}
    }
    for $i (keys %{ $PARMS }) {
	delete $PARMS->{$i} if ($i !~ /ATOMTYPES|VDW/);
    }
	
}    

sub createLammpsFiles {
    my ($data,  $fflist) = @_;
    my ($curr, $ffs, $execStr, $bgfFile, $saveName);

    system("mkdir -p _lmptmp");
    $execStr = "/home/yjn1818/scripts/createLammpsInput.pl ";
    for $curr (@{ $fflist }) {
	$saveName = abs_path($curr);
	$ffs .= "$saveName ";
    }

    for $curr (keys %{ $data }) {
	$bgfFile = abs_path($data->{$curr}{FILE});
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$//;
	chdir "_lmptmp";
	if (system("${execStr} -f \"$ffs\" -b $bgfFile -s $saveName> /dev/null")) {
	    die "ERROR: Cannnot execute \"${execStr} -f \"$ffs\" -b $bgfFile\"\n";
	}
	$data->{$curr}{LMPDAT} = "_lmptmp/data.${saveName}";
	chdir "../";
    }
}

sub getOptimalParams {
    my ($parms, $outfile) = @_;
    my ($tol, $iter, $opt, $error, $i, $j);

    $tol = 5e-02;
    $iter = 1000;

    $outfile =~ s/\.\w+$//;
    $outfile .= "_residuals.dat";
    open $DATFILE, "> $outfile" or die "ERROR: Cannot write to $outfile: $!\n";
    #### SIMPLEX DOWNHILL ALGORITHM #####
    ($opt, $error) = MinimiseND($parms->{OPT}{guess}, $parms->{OPT}{scale}, \&calcResiduals, $tol, $iter);
    print "Final Params\n";
    for $i (0 .. $#{ $opt }) {
	printf "%-10s %2s: %.5f\n", $parms->{OPT}{map}[$i]{label}, $parms->{OPT}{map}[$i]{parms}{name}, $opt->[$i];
    }
    print $DATFILE "#ERROR: $error\n";
    close $DATFILE;

}

sub calcResiduals {
    my ($execStr, $i, $pair_coeff, $inp, $datfile, $valid, $tot);

    $count++;
    $valid = $tot = 0;

    $inp = "_lmptmp/in.tmp";
    system("cp $lammpsInput $inp");
    $execStr = "$lammpsBinary -in $inp -var datafile ";
    $pair_coeff = getPairCoeffs($PARMS, \@_);
    open DAT, ">> $inp" or die "ERROR: Cannot write to $inp: $!\n";
    print DAT $pair_coeff;
    close DAT;

    for $i (keys %{ $DATA }) {
        $datfile = $DATA->{$i}{LMPDAT};
        open LMP, "${execStr} $datfile |" or die "Cannot execute \"${execStr} $datfile\": $!\n";
        $valid = 0;
        while (<LMP>) {
            chomp;
            if ($_ =~ /interact\s+=\s+(\-?\d+\.\d+)/i) {
                $tot += sprintf("%.5f", ($DATA->{$i}{WEIGHT}**2 * ($1 - $DATA->{$i}{ENG})**2));
                $valid = 1;
            } elsif ($_ =~ /forces\s+=\w+(\d+\.\d+)/) {
		$tot += sprintf("%.5f", ($DATA->{$i}{WEIGHT}**2 * $1));
		$valid = 1;
	    }
        }
        close LMP;
        die "ERROR: No valid data found file executing \"$execStr $datfile\"\n" if (! $valid);
    }

    print $DATFILE "$count $tot\n";
    printf "%10d%12.5f\n",$count,$tot;
    return $tot;
}

sub getPairCoeffs {
    my ($ffData, $newParms) = @_;
    my ($outStr, $i, $j, @vdws, @tmp, $curr); 
    my ($pair, $vdwType, $type1, $type2);
    
    if ($newParms) {
	for $i (0 .. $#{ $ffData->{OPT}{map} }) {
	    next if (! defined($newParms->[$i]));
	    $j = $ffData->{OPT}{map}[$i]{name}; # j points to $PARMS{VDW}{type1}{type2}
	    $j->{1}{VALS}[$ffData->{OPT}{map}[$i]{parms}{val}] = $newParms->[$i]; # update j val to reflect new val from simplex
	}
    }

    for $i (keys %{ $ffData->{VDW} }) {
	for $j (keys %{ $ffData->{VDW}{$i} }) {
	    $type1 = $ffData->{ATOMTYPES}{$i}{INDEX};
	    $type2 = $ffData->{ATOMTYPES}{$j}{INDEX};
	    ($type1, $type2) = ($type2, $type1) if ($type1 > $type2);
	    $curr = $ffData->{VDW}{$i}{$j}{1};
	    $vdwType = uc $curr->{TYPE};
	    @vdws = @{ $curr->{VALS} };
	    if ($vdwType =~ /^MORSE|STRETCH_MORSE$/) {
		($vdws[0], $vdws[1]) = ($vdws[1], $vdws[0]);
		($vdws[1], $vdws[2]) = ($vdws[2], $vdws[1]);
		$vdws[1] /= ($vdws[2] * 2); # changed from div to multi 07/28/2007
		$vdws[3] /= ($vdws[2] * 2) if ($vdwType =~ /^STRETCH_MORSE/); #aplha2 for stretch morse
	    } elsif ($vdwType =~ /^LJ_6_12$/) {
		($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
		$vdws[1] = $vdws[1] / (2**(1/6));
		if ($#vdws > 1) {
		    ($vdws[2],$vdws[3]) = ($vdws[3],$vdws[2]);
		    $vdws[3] = $vdws[3] /(2**(1/6));
		}
	    } elsif ($vdwType =~ /^EXPO_6$/) {
		($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
		@tmp = @vdws;
		@vdws = ();
		$vdws[0] = $tmp[0] * (6/($tmp[2]-6)) * exp($tmp[2]);
		$vdws[1] = $tmp[1]/$tmp[2];
		$vdws[2] = $tmp[1]**6 * $tmp[0] * ($tmp[2]/($tmp[2] - 6));
	    } elsif ($vdwType =~ /^BUCKINGHAM/) {
		($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
		$vdws[1] = 1/$vdws[1];
	    }
	    $pair->{$type1}{$type2}{TYPE} = $curr->{TYPE};
	    @{ $pair->{$type1}{$type2}{VALS} } = @vdws;
	}
    }

    for $type1 (keys %{ $pair }) {
	for $type2 (keys %{ $pair }) {
	    next if ($type1 > $type2);
	    $pair->{$type1}{$type2}{VALS} = mix($type1, $type2, $pair) if (! exists($pair->{$type1}{$type2}));
	    $outStr .= sprintf("%-12s%8d%8d    @{ $pair->{$type1}{$type2}{VALS} }\n", "pair_coeff", $type1, $type2);
	}
    }
    $outStr .= "\ncompute         interact group1 group/group group2\n";
    #$outStr .= "compute         forceX all reduce sum fx\n";
    #$outStr .= "compute         forceY all reduce sum fy\n";
    #$outStr .= "compute         forceZ all reduce sum fz\n";
    #$outStr .= "variable        forces equal c_forceX^2+c_forceY^2+c_forceZ^2\n";
    $outStr .= "thermo_style    custom step cpu pe c_interact evdwl\n";
    $outStr .= "thermo_modify   line multi\n";
    $outStr .= "run             0\n";

    return $outStr;
}

sub mix  {
    my ($t1, $t2, $pairData) = @_;
    my ($aa, $bb, @vals);

    $aa = $pairData->{$t1}{$t1};
    $bb = $pairData->{$t2}{$t2};
    return () if ($aa->{TYPE} ne $bb->{TYPE});
    
    if (uc($aa->{TYPE}) eq "VDW_MORSE") {
	$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
	$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
	$vals[2] = 0.5*($aa->{VALS}[2]+$bb->{VALS}[2]);
    } elsif (uc($aa->{TYPE}) eq "LJ_6_12") {
	$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
	$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
	if ($#{ $aa->{VALS} } > 1) {
	    $vals[2] = sqrt($aa->{VALS}[2]*$bb->{VALS}[2]);
	    $vals[3] = 0.5*($aa->{VALS}[3]+$bb->{VALS}[3]);
	}
    } elsif (uc($aa->{TYPE}) =~ /EXPO_6|BUCKINGHAM/) {
	$vals[0] = sqrt($aa->{VALS}[0]*$bb->{VALS}[0]);
	$vals[1] = 0.5*($aa->{VALS}[1]+$bb->{VALS}[1]);
	$vals[2] = sqrt($aa->{VALS}[2]*$bb->{VALS}[2]);
    }
    return \@vals;
}


sub assignBGFS {
    my ($data, $bgfs) = @_;
    my ($i, $isValid, $j, $count);

    $count = 0;
    for $i (keys %{ $data }) {
	$isValid = 0;
	for $j (0 .. $#{ $bgfs }) {
	    if ($bgfs->[$j] =~ /$i/) {
		$data->{$i}{FILE} = $bgfs->[$j];
		splice @{ $bgfs }, $j, 1;
		$isValid = 1;
		last;
	    }
	}
	die "ERROR: No corresponding bgf file found while search for $i!\n"
	    if (! $isValid);
	$count++;
    }
    print "read in $count datavalues and assigned corresponding bgffiles...";
}

sub getObjectiveFunctionData {
    my ($eFile) = $_[0];
    my (%EFUNC, $rec);

    open ENGFILE, $eFile or die "ERROR: Cannot open energy file $eFile: $!\n";
    while (<ENGFILE>) {
	chomp;
	if ($_ =~ /^\s*(\-?\d+\.?\d*)\s+(\-?\d+\.\d+)\s+(\d+\.?\d*E?\-?\d*)/) {
	    $rec = (
		    {
			"DISP"   => $1,
			"ENG"    => $2,
			"WEIGHT" => $3,
		    }
		    );
	    $EFUNC{$1} = $rec;
	} elsif ($_ =~ /^Constrain\s+(\w+)\s+(\w+)\s+(\w+)/) {
	    $FREEZE->{$1}{$2}{$3} = 1;
	    $FREEZE->{$2}{$1}{$3} = 1;
	}
    }
    close ENGFILE;
    die "ERROR: Energy file $eFile does not contain any valid information!\n"
	if (! %EFUNC);

    return \%EFUNC;
}

sub getBounds {
    my ($vdwData, $constraints, $optScale) = @_;
    my ($fMap, $pGuess, $pScale, $pOrder, $i, $j, $k, $count, $factor);

    $count = 0;
    $factor = 0.05;
    for $i (keys %{ $vdwData }) {
	for $j (keys %{ $vdwData->{$i} }) {
	    next if ($i lt $j);
	    $pOrder = getParmOrder($vdwData->{$i}{$j}, $optScale);
	    for $k (@{ $pOrder->{parms} }) {
		next if (exists($constraints->{$i}{$j}{$k->{name}}));
		push @{ $pGuess }, $vdwData->{$i}{$j}{1}{VALS}[$k->{val}];
		push @{ $pScale }, $factor * $vdwData->{$i}{$j}{1}{VALS}[$k->{val}];
		push @{ $fMap }, {"name" => $pOrder->{name}, "parms" => $k, "label" => "${i}-${j}" };
	    }
	}
    }

    die "ERROR: No valid info found while creating objective function!\n" if (! $fMap);

    return ($pGuess, $pScale, $fMap);
}

sub getParmOrder {
    my ($parm, $optScale) = @_;
    my ($rec);

    $rec->{name} = $parm;
    $rec->{parms}[0] = {"name" => "r0", "val" => 0,};
    $rec->{parms}[1] = {"name" => "d0", "val" => 1,};
    push @{ $rec->{parms} }, ({"name" => "s1", "val" => 2}) if ($optScale && $#{ $parm->{1}{VALS} }> 1); #morse/exp6
    push @{ $rec->{parms} }, ({"name" => "s2", "val" => 3}) if ($optScale && $#{ $parm->{1}{VALS} }> 2); #stretched morse/exp6
  
    return $rec;
}

sub init {
    my (%OPTS);
    my ($tmp, $findCmd, $ffList, $singleBGF, $flist, $rec);

    $createLammpsScript = "/home/yjn1818/scripts/createLammpsInput.pl";
    $optScale = 0;
    getopt('dflsiao',\%OPTS);
    for ("d", "f", "i") {
	die "usage: $0 -d \"bgf file(s):energy file ...\" -f forcefield -i lammps input file\n\t\t" . 
	    "-l (lammps binary) -s (save name) -a (optimize vdw scale = no) -o (opt off diag = no)\n"
	    if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($flist, $ffList, $lammpsBinary, $saveName, $lammpsInput, $optOffDiag) = 
	($OPTS{d}, $OPTS{f}, $OPTS{l}, $OPTS{s}, $OPTS{i}, $OPTS{o});
    $lammpsBinary = "/home/yjn1818/programs/bin/lmp_serial64" if (! defined($lammpsBinary));
    #die "ERROR: Cannot locate lammps binary $lammpsBinary: $!\n" if (! -e $lammpsBinary);
    $optScale = 1 if (exists($OPTS{a}) and $OPTS{a} =~ /1|yes/i);
    while ($flist =~ /(\S+):(\S+)/g) {
	$rec = { "ENG" => $2, "BGF" => $1, };
	push @{ $FILES }, $rec;
    }
    die "ERROR: Expected \"bgf file(s):energy file ...\"n for -d, got \"$flist\"\n" if (! defined($FILES));
    FileTester($lammpsInput);
    while ($ffList =~ /(\S+)/g) {
	push @{ $ffData }, $1 if (-e $1 and ! -d $1 and -r $1 and -T $1);
    }
    die "ERROR: No valid ffs found while searching \"$ffList\"\n" if (! $ffData);
    if (! defined($saveName)) {
	$saveName = basename($ffData->[0]);
	$saveName =~ s/\.\w+$//;
	$saveName .= "_opt.ff";
    }
    $prefix = $saveName;
    $prefix =~ s/\.\w+$//;
    if (defined($optOffDiag)) {
	undef($optOffDiag) if ($optOffDiag !~ /^1|0$/);
    }
    print "Done\n";
}
