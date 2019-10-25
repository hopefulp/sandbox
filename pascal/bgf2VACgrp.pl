#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use Packages::General qw(FileTester);
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::ManipAtoms qw(GetMols);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub createGrpFile;
sub getGrps;
sub numerically { ($a<=>$b); }
sub getConstraints;

my ($field, $bgfFile, $saveName, $shake, $opts);
my ($ATOMS, $BONDS, $GRPS, $zeroConstraints);

$|++;
&init;

print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 1);
&GetMols($ATOMS, $BONDS) if ($field =~ /MOL/i);
print "Done\n";
$GRPS = getGrps($ATOMS, $field);
print "Creating Group file ${saveName}...";
&createGrpFile($GRPS, $saveName, $shake, $opts);
print "Done\n";

sub getGrps {
    my ($atoms, $field) = @_;
    my ($gHASH, $i, $hKey, $tmp, $gMap, @list);

    for $i (keys %{ $atoms }) {
	$tmp->{ $atoms->{$i}{$field} } = 1 if ($field !~ /MOL/i);
	$tmp->{ ${ $atoms->{$i}{$field} } } = 1 if ($field =~ /MOL/i);
    }
    if ($field =~ /CHAIN|FFTYPE|RESNAME/) {
	@list = sort {$a cmp $b} keys %{ $tmp };
    } else {
	@list = sort numerically keys %{ $tmp };
    }
    for $i (0 .. $#list) {
	$gMap->{$list[$i]} = $i;
    }
    for $i (keys %{ $atoms }) {
	$hKey = $gMap->{ $atoms->{$i}{$field} } if ($field !~ /MOL/i);
	$hKey = $gMap->{ ${ $atoms->{$i}{$field}  }} if ($field =~ /MOL/i);
	$gHASH->{$hKey}{$i} = 1;
    }

    return $gHASH;
}


sub createGrpFile {
    my ($data, $saveFile, $shakeOpt, $otherOpts) = @_;
    my ($i, @tmp, @tmp2, $j, $constraints);
    my ($groups, $start, $prev, $counter);

    @tmp = sort numerically keys %{ $data };
    $i = scalar(@tmp);
    open OUTDATA, "> $saveFile" || die "ERROR: Cannot write to ${saveFile}: $!\n";
    print OUTDATA "Total Groups: $i\n";

    for $i (@tmp) {
	@tmp2 = sort numerically keys %{ $data->{$i} };
	$j = scalar(@tmp2);
	print OUTDATA "Group " . ($i+1) . " Atoms $j\n";
        $start = $prev = -1;
        $groups = "";
        $counter = 0;
        for $j (@tmp2) {
            if ($start == -1) {
                $start = $j;
            } elsif (($j - $prev) > 1) {
		if (($prev - $start) > 1) {
		    $groups .= "${start} - ${prev} ";
		} elsif (($prev - $start) > 0) {
		    $groups .= "${start} ${prev} ";
		} else {
		    $groups .= "${start} ";
		}
                $start = $j;
                $counter++;
            }
            $prev = $j;
            if ($counter == 5) {
                $counter = 0;
                $groups = substr($groups, 0, -1);
                $groups .= "\n";
            }
        }
        $groups .= "${start} - ${prev} ";
        $groups = substr($groups, 0, -1);
        print OUTDATA "$groups\n";
    }
    if ($shakeOpt) { 
	$constraints = getConstraints($data, $otherOpts, $zeroConstraints);
    } else { 
	$constraints = "";
    }
    print OUTDATA $constraints;
    close OUTDATA;
}

sub getConstraints {
    my ($data, $oOpts, $setZero) = @_;
    my ($i, $j, $mols, @tmp, $cStr, $rStr, $lStr, $count, $eType, $atom, $tot);
    
    @tmp = sort numerically keys %{ $data };

    $cStr = "Constraints\n";
    $rStr = "RotationalSymmetryNumber\n";
    $lStr = "LinearMoleculeFlag\n";
    for $i (@tmp) {
	$count = 0;
	$tot = 0;
	$mols = ();
	$mols = GetMols($ATOMS, $BONDS, $data->{$i});
	for $j (keys %{ $mols }) {
	    for $atom (keys %{ $mols->{$j}{MEMBERS} }) {
		$count++ if ($ATOMS->{$atom}{FFTYPE} =~ /^H/);
	    }
	    $tot += $mols->{$j}{MOLSIZE};
	    $count++ if ($mols->{$j}{MOLSIZE} == 3);
	}
	$cStr .= "$count " if (! $setZero);
        $cStr .= "0 " if ($setZero);
	if ((scalar(keys %{ $mols }) > 0) and (($tot/scalar(keys %{ $mols })) == 3)) {
	    $rStr .= "2 "; 
	} else {
	    $rStr .= "1 ";
	}
	$lStr .= "0 ";
    }
    $rStr = $lStr = "" if(! $oOpts);
    return "$cStr\n$rStr\n$lStr\n";
}

sub init {
    my (%OPTS);

    getopt('fbscoz',\%OPTS);

    ($field, $bgfFile, $saveName, $shake, $opts, $zeroConstraints) = ($OPTS{f},$OPTS{b},$OPTS{s}, $OPTS{c}, $OPTS{o}, $OPTS{z});

    die "usage:$0 -b bgf file -f (field=chain|resid|resname|mol|fftype) -c (has shake=no) -o (other opts=yes if shake) -z (zero constraints) -s (savename)\n" 
	if (! defined($OPTS{b}));
    print "Initializing...";

    FileTester($bgfFile);
    $field = "CHAIN" if (! defined($field));
    if ($field =~ /(chain|resname|res|molsize|mol|fftype)/i) {
	$field = uc $1;
	$field .= "ECULEID" if ($field eq "MOL");
	$field .= "NUM" if ($field eq "RES");
    }

    if (! defined($saveName)) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$//;
	$saveName .= ".grps";
    }
    $shake = 0 if (! defined($shake) or $shake !~ /1|yes/i);
    $shake = 1 if (defined($shake) and $shake =~ /1|yes/i);
    $opts = 1 if (! defined($opts) or $opts !~ /0|no/i);
    $opts = 0 if (defined($opts) and $opts =~ /0|no/i);
    $zeroConstraints = 0 if (! defined($zeroConstraints) or $zeroConstraints !~ /1|yes/i);
    $zeroConstraints = 1 if ($zeroConstraints =~ /1|yes/i);
    print "Done\n";
}
