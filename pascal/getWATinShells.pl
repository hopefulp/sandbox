#!/usr/bin/perl
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use warnings;
no warnings "recursion";
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use Packages::General qw(FileTester IsInteger CoM);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use Packages::ManipAtoms qw(SplitAtomsByMol ImageAtoms);
use Packages::BOX qw(GetBox);

sub init;
sub numerically { ($a<=>$b); };
sub findWatShells;
sub createFile;
sub showUsage;
sub getShakeDOF;
sub getMolTypes;
sub updateBGF;
sub storeCOORDS;
sub findCloseAtoms;
sub reImage;

my ($bgfFile, $cerius2FF, $shellStr, $saveName, $watShells);
my ($ATOMS, $BONDS, $HEADERS, $MOLS, $MOLtypes, $SHELLS, $num, $BOX);

$|++;
$watShells = &init;
print "Getting atom data from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$MOLS = SplitAtomsByMol($ATOMS);
$BOX = GetBox($ATOMS, undef, $HEADERS);
print "Done\nDetermining types of molecules...";
($MOLtypes, $num) = getMolTypes($ATOMS, $MOLS);
&reImage($ATOMS, $MOLtypes, $BOX);
print "Done\nGetting water in shells...";
my $start = time();
$SHELLS = findWatShells($MOLtypes,$watShells, $ATOMS, $num);
my $end = time();
print "Done\nCreating $saveName...";
&createFile($ATOMS, $SHELLS, $saveName);
$saveName = &updateBGF($ATOMS, $saveName);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\nElapsed " . ($end-$start). " s\n";

sub updateBGF {
    my ($atoms, $savePrefix) = @_;
    my ($i, $j);

    $savePrefix =~ s/\.\w+$//;
    $savePrefix .= ".bgf";
    for $i (keys %{ $atoms }) {
	next if (! exists($atoms->{$i}{SHELL}));
	$atoms->{$i}{CHAIN} = chr(64 + $atoms->{$i}{SHELL});
     }

     return $savePrefix;
}

sub createFile {
    my ($atoms, $wshells, $save) = @_;
    my (@tmp, @tmp2, $i, $j, $start, $prev, $groups, $counter, $dofStr);
    
    @tmp = keys %{ $wshells };
    open OUTDATA, "> $save" || die "ERROR: Cannot create $save: $!\n";
    print OUTDATA "Total Groups: " . ($#tmp + 1) . "\n";

    @tmp = sort numerically keys %{ $wshells };
    $counter = 0;
    for $i (@tmp) {
	@tmp2 = sort numerically keys %{ $wshells->{$i}{ATOMS} };
        $dofStr .= getShakeDOF($wshells->{$i}{ATOMS}, $atoms) . " ";
	print OUTDATA "Group $i Atoms " . $wshells->{$i}{COUNT} . "\n";
	$start = $prev = -1;
	$groups = "";
	$counter = 0;
	for $j (@tmp2) {
	    if ($start == -1) {
		$start = $j;
	    } elsif (($j - $prev) > 1) {
		$groups .= "${start} - ${prev} ";
		$start = $j;
		$counter++;
	    }
	    $prev = $j;
	    if ($counter == 10) {
		$counter = 0;
		$groups = substr($groups, 0, -1);
		$groups .= "\n";
	    }
	}
	$groups .= "${start} - ${prev} ";
	$groups = substr($groups, 0, -1);
	print OUTDATA "$groups\n";
    }
    print OUTDATA "Constraints\n$dofStr\n";
    close OUTDATA;
}

sub findWatShells {
    my ($system, $shellData, $atomData, $offset) = @_;
    my ($mol, $DIST, $dSquared, $solvAtom, $soluAtom, $tot); 
    my (@dSort, $dist, $i, $shell, $SHELLS, $CoM, @soluList);

    $tot = scalar(@{ $shellData }) + $offset;
    for $mol (keys %{ $system->{SOLV}{LIST} }) {
	$DIST = ();
	$CoM = ();
	$i = 0;
	for $solvAtom (keys %{ $system->{SOLV}{LIST}{$mol} }) {
	    $CoM->{XCOORD} += $atomData->{$solvAtom}{XCOORD};
	    $CoM->{YCOORD} += $atomData->{$solvAtom}{YCOORD};
	    $CoM->{ZCOORD} += $atomData->{$solvAtom}{ZCOORD};
	    $i++;
	}
        $CoM->{XCOORD} /= $i;
	$CoM->{YCOORD} /= $i;
	$CoM->{ZCOORD} /= $i;
	@soluList = keys %{ $system->{SOLU}{LIST} };
	for $soluAtom (@soluList) {
	    $dist = ($atomData->{$soluAtom}{XCOORD} - $CoM->{XCOORD})**2 +
	    	    ($atomData->{$soluAtom}{YCOORD} - $CoM->{YCOORD})**2 +
		    ($atomData->{$soluAtom}{ZCOORD} - $CoM->{ZCOORD})**2;
	    $DIST->{$dist} = $soluAtom;
	}
	undef $shell;
	@dSort = sort numerically keys %{ $DIST };
	$dist = sqrt($dSort[0]);
	for $i (0 .. $#{ $shellData }) {
	    if ($shellData->[$i] > $dist) {
		$shell = ($i + $offset);
		last;
	    }
	}
	$shell = $tot if (! defined($shell));
        for $solvAtom (keys %{ $system->{SOLV}{LIST}{$mol} }) {
	    $SHELLS->{$shell}{ATOMS}{$solvAtom} = 1;
	    $SHELLS->{$shell}{COUNT}++;
	    $atomData->{$solvAtom}{SHELL} = $shell;
	}
    }

    $SHELLS->{1}{ATOMS} = $system->{SOLU}{LIST};
    $SHELLS->{1}{COUNT} = $system->{SOLU}{NUM};
    if ($offset == 3) {
	$SHELLS->{2}{ATOMS} = $system->{IONS}{LIST};
	$SHELLS->{2}{COUNT} = $system->{IONS}{NUM};
    }

    return $SHELLS;
}


sub getMolTypes {
    my ($atoms, $mols) = @_;
    my ($i, $j, $count, $molSize, $molStats, @tmp, $atom, $types, $offset);

    $offset = 2;
    for $i (keys %{ $mols } ) {
	@tmp = keys %{ $mols->{$i} };
	$atom = $tmp[0];
	$molSize = $atoms->{$atom}{MOLECULE}{SIZE};
	$molStats->{$molSize}{LIST}{$i} = $mols->{$i};
	$molStats->{$molSize}{NUM}++;
    }

    if (exists($molStats->{1})) { #ions
	$types->{IONS}{NUM} = $molStats->{1}{NUM};
        for $i (keys %{ $molStats->{1}{LIST} }) {
	    for $j (keys %{ $molStats->{1}{LIST}{$i} }) {
		$types->{IONS}{LIST}{$j} = 1;
	    }
	}
	print $molStats->{1}{NUM} . " ions...";
	delete $molStats->{1};
	$offset++;
    }
    if (exists($molStats->{3})) { #waters
	$types->{SOLV} = $molStats->{3};
	print $molStats->{3}{NUM} . " waters...";
	delete $molStats->{3};
    } else {
	die "ERROR: No water atoms found!\n";
    }

    $count = 0;
    for $i (keys %{ $molStats }) {
	for $j (keys %{ $molStats->{$i}{LIST} }) {
	    for (keys %{ $molStats->{$i}{LIST}{$j} }) {
		$types->{SOLU}{LIST}{$_} = 1;
		$count++;
	    }
	}
    }
    $types->{SOLU}{NUM} = $count;
    die "ERROR: No solute atoms found!\n" if (! $count);
    print "$count solute atoms...";

    return ($types, $offset);
}


sub init {
    my (@tmp, $i, %SHELL, $counter, %OPTS, $usage);

    getopt('bsw',\%OPTS);

    $usage = &showUsage;
    for ("b", "w") {
	die "$usage\n" if (! $OPTS{$_});
    }
    
    print "Initializing...";
    ($bgfFile, $shellStr, $saveName) = ($OPTS{b},$OPTS{w},$OPTS{s});
    FileTester($bgfFile);
    
    if ($shellStr =~ /\s+/) {
	@tmp = split /\s+/, $shellStr;
    } else {
	$tmp[0] = $shellStr;
    }
    $counter = 0;
    for $i (@tmp) {
	if ($i =~ /(\d+\.*\d+)/) {
	    $SHELL{$counter} = $1;
	    $counter++;
	}
    }
    die "Expected at least one decimal(integer) value for the shell limts. Got \"$shellStr\"\n" 
	if (! $counter);
    
    if (! defined($saveName)) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$/_watshell\.dat/;
    }
    print "Done\n";
    return \@tmp;
}

sub reImage {
    my ($atoms, $molTypes, $box) = @_;
    my ($molData, $i, $j, $center, $molList, @tmp);

    #move solute to center of box
    for $i (keys %{ $molTypes->{SOLU}{LIST} }) {
	$molData->{$i} = $atoms->{$i};
    }
    $center = CoM($molData); 
    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    $box->{XCOORD} = $box->{X};
    $box->{YCOORD} = $box->{Y};
    $box->{ZCOORD} = $box->{Z};
    for $j (@tmp) {
        $box->{$j}{CENTER} = $box->{$j}{len}/2;
        for $i (keys %{ $atoms }) {
           $atoms->{$i}{$j} += ($box->{$j}{CENTER} - $center->{$j});
        }
    }

    #reimage ions
    for $i (keys %{ $molTypes->{IONS}{LIST} }) {
	$molData = ();
	$molData->{$i} = $atoms->{$i};
	$center = $atoms->{$i};
	ImageAtoms($molData, $center, $box);
    }

    #reimage solvent/water
    for $i (keys %{ $molTypes->{SOLV}{LIST} }) {
	$molData = ();
	for $j (keys %{ $molTypes->{SOLV}{LIST}{$i} }) {
	    $molData->{$j} = $atoms->{$j};
	}

        $center = CoM($molData);
	ImageAtoms($molData, $center, $box);
    }
}


sub getShakeDOF {
    my ($select, $atoms) = @_;
    my ($i, $dof);
    $dof = 0;
    for $i (keys %{ $select }) {
        if ($atoms->{$i}{FFTYPE} =~ /^H/) { #if this is a hydrogen
            $dof++;
            if ($atoms->{$i}{MOLSIZE} == 3) {
                $dof += 0.5;
            }
        }
    }
    return $dof;
}

sub showUsage {
    my ($usage) = "This script will break up the bgf file into:  solute, ions, water shell 1...n, and bulk\n\n" . 
	"usage: $0 -b bgf file -w water shell distance(s) -s [save name]\n" .
	"\toptions:\n" .
	"\t-b bgf file: the location of the bgf file\n" .
	"\t-w water shell distances(s): the distance from the surface of the 1st water shell.\n" .
	"\t    you can specify multiple water shell distances by enclosing them in \"\" qoutes\n" .
	"\t-s [save name]: (optional) the name to save the group information generated by the script\n" .
	"\t    if not specified, will be bgf_file.grps. Each water shell will be save_name_shellx.bgf\n";
    return $usage;
}
