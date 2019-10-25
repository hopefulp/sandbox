#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF createHeaders addBoxToHeader);
use Packages::General qw(FileTester CombineMols CoM);
use Packages::BOX qw(GetBox);
use Packages::ManipAtoms qw(GetMols CenterSystem);

sub init;
sub replicateCell;
sub transMol;
sub getPBCbonds;
sub updatePBCbonds;
sub updateBond;

my ($replicate, $bgfFile, $saveFile, $doCenter, $bondImages);
my ($CENTER, $ATOMS, $BONDS, $HEADERS, $BOX);
$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$BOX = GetBox($ATOMS, undef, $HEADERS);
print "Done\nReplicating cell $replicate->{STRING}...";
($ATOMS, $BONDS, $BOX) = replicateCell($ATOMS, $BONDS, $BOX, $replicate, $doCenter, $bondImages);
&GetMols($ATOMS, $BONDS);
print "Done\nCreating BGF file $saveFile...";
$HEADERS = createHeaders(undef, $saveFile);
#$CENTER = CoM($ATOMS);
#for (keys %{ $BOX }) {
    #$CENTER->{$_ . "COORD"} -= (($BOX->{$_}{hi} + $BOX->{$_}{lo})/2 + $BOX->{$_}{lo} - $CENTER->{$_ . "COORD"});
#}
#&CenterSystem($ATOMS, $CENTER, ());
&addBoxToHeader($HEADERS, $BOX);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub transmol {
    my ($unit, $box, $dim, $disp) = @_;
    my ($i,%ATOMS, $offset);
    
    $offset = $box->{$dim}{len} * $disp;
    for $i (keys %{ $unit } ) {
	%{ $ATOMS{$i} } = %{ $unit->{$i} };
    }
    for $i (keys %{ $unit }) {
	$ATOMS{$i}{$dim . "COORD"} += $offset;
    }
    return \%ATOMS;
}

sub replicateCell {
    my ($atoms, $bonds, $box, $dims, $centerMol, $createNewbonds) = @_;
    my ($i, $j, $cellAtoms, $cellBonds, $unitAtoms, $pbcbonds, $offset, $tot);
    
    if ($centerMol) {
	for $i ("X", "Y", "Z") {
	    $j = -1 * int(($dims->{$i} -1)/2);
	    for ($j .. -1) {
		$atoms = transmol($atoms, $box, $i, $_);
	    }
	}
    }

    for $i ("X", "Y", "Z") {
	$unitAtoms = ();
	$cellBonds = ();
	$pbcbonds = ();
	$pbcbonds = getPBCbonds($atoms, $bonds, $box) if ($createNewbonds);
	$tot = 0;
	$offset = 0;
	for $j (keys %{ $atoms }) {
	    %{ $unitAtoms->{$j} } = %{ $atoms->{$j} };
	    @{ $cellBonds->{$j} } = @{ $bonds->{$j} };
	    $tot++;
	}
	for $j (1 .. ($dims->{$i} - 1)) {
	    $offset += $tot;
	    $cellAtoms = transmol($unitAtoms, $box, $i, $j);
	    ($atoms, $bonds) = CombineMols($atoms, $cellAtoms, $bonds, $cellBonds);
	    updatePBCbonds($atoms, $bonds, $pbcbonds->{$i}, $offset, $tot, $i) if(exists($pbcbonds->{$i}));
	}
	$box->{$i}{hi} = ($box->{$i}{len} * $dims->{$i});
	$box->{$i}{lo} = 0;
	$box->{$i}{len} = ($box->{$i}{len} * $dims->{$i});
    }
    return ($atoms, $bonds, $box);
}

sub updatePBCbonds {
    my ($atoms, $bonds, $bondlist, $atomOffset, $tot, $dim) = @_;
    my ($i, $j, $atom1, $atom2, $atom3, $atom4, @list, $pos);

    for $i (keys %{ $bondlist }) {
	@list = keys %{ $bondlist->{$i} };
	for $j (0 .. $#list) {
	    $pos = $bondlist->{$i}{ $list[$j] };
	    $atom1 = $i;
	    $atom2 = $list[$j];
	    $atom3 = $atom1 + $atomOffset;
	    $atom4 = $atom2 + $tot;
	    #for atom1 , delete bond to atom2 and form bond to atom4
	    $bonds->{$atom1}[$pos->[0]] = $atom4;
            #for atom2, delete bond to atom1 and form bond to atom3
	    $bonds->{$atom2}[$pos->[1]] = $atom3;
	    delete $atoms->{$atom2}{"DISP${dim}"};
	    #$atoms->{$atom2}{"DISP${dim}"}[$pos->[1]] = 0;
            #for atom3, delete bond to atom4 and form bond to atom2
            $bonds->{$atom3}[$pos->[0]] = $atom2;
	    delete $atoms->{$atom3}{"DISP${dim}"};
            #$atoms->{$atom3}{"DISP${dim}"}[$pos->[0]] = 0;
            #updateBond($atoms->{$atom3}, $bonds->{$atom3}, $atom4, $atom2);
            #delete $atoms->{$atom3}{"DISP${dim}"};
            #for atom4, delete bond to atom3 and form bond to atom1
            $bonds->{$atom4}[$pos->[1]] = $atom1;
            #updateBond($atoms->{$atom4}, $bonds->{$atom4}, $atom3, $atom1);
            #now update bondlist
            delete $bondlist->{$i}{$atom2};
            $bondlist->{$i}{$atom4}[0] = $pos->[0];
	    $bondlist->{$i}{$atom4}[1] = $pos->[1];
	}
    }
    print "";
}

sub updateBond {
    my ($atomData, $atomBonds, $bondAtom, $newAtom) = @_;
    my ($i);

    $i = 0;
    while($i<= $#{ $atomBonds }) {
	if ($atomBonds->[$i] == $bondAtom) {
	     $atomBonds->[$i] = $newAtom;
	     last;
	} else {
	    $i++;
	}
    }
}

sub getPBCbonds {
    my ($atoms, $bonds, $box) = @_;
    my ($i, $j, $k, $l, $dist, $atom1, $atom2, $bondlist, $sign,$flag2);

    #search for multiple bond entries to same atom and image flags
    for $i (keys %{ $bonds }) {
	$atom1 = $i;
        for $j (0 .. $#{ $bonds->{$i} }) {
            for $k ("X", "Y", "Z") {
		$sign = 0;
		$sign = $atoms->{$i}{"DISP${k}"}[$j] if (exists($atoms->{$i}{"DISP${k}"}));
		next if($sign>-1);
		#find symmetric entry
		undef($flag2);
		$atom2 =  $bonds->{$i}[$j];
		for $l (0 .. $#{ $bonds->{$atom2} }) {
		    if ($bonds->{$atom2}[$l] == $atom1) {
			if (exists($atoms->{$atom2}{"DISP${k}"}) and $atoms->{$atom2}{"DISP${k}"}[$l] == -1*$sign) {
			    $flag2 = $l;
			    last;
			}
		    }
		}
		next if (! defined($flag2));
		$bondlist->{$k}{$i}{$bonds->{$i}[$j]} = ([$j, $flag2]);
            }
        }
    }

    return $bondlist if(defined($bondlist));
    #search for bonds by distance
    for $i (keys %{ $bonds }) {
        $atom1 = \%{ $atoms->{$i} };
	for $j (0 .. $#{ $bonds->{$i} }) {
            $atom2 = \%{ $atoms->{$bonds->{$i}[$j]} };
	    for $k ("X","Y","Z") {
		$dist = $atom1->{"${k}COORD"} - $atom2->{"${k}COORD"};
		next if ($dist == 0);
		$sign = $dist/abs($dist);
		if ($sign < 0) {
		    $dist *= -1;
		}
		if ($dist > 5 and $dist > $box->{$k}{len}/3) { #search for long bonds
		    $bondlist->{$k}{$bonds->{$i}[$j]}{$i} = ([$j,$sign]);
		}
	    }
	}
    }

    return $bondlist;
}

sub init {
    my (%OPTS, $rStr);
    getopt('bdsci',\%OPTS);
    die "usage: $0 -b bgf file -d \"x y z\" replication\n" . 
	"\t-s (save name) -c (center=no) -i (make periodic bonds=yes)\n" 
	if (! exists($OPTS{b}) or ! exists($OPTS{d}));
    print "Initializing...";
    ($bgfFile, $rStr, $saveFile, $doCenter, $bondImages) = ($OPTS{b}, $OPTS{d}, $OPTS{s}, $OPTS{c}, $OPTS{i});
    FileTester($bgfFile);
    if ($rStr !~ /(\d+)\s+(\d+)\s+(\d+)/) {
	die "ERROR: Expected integers for x,y and z replication. Got \"$rStr\"\n";
    } else {
	$replicate = (
		      {
			  "X"      => $1,
			  "Y"      => $2,
			  "Z"      => $3,
			  "STRING" => "${1}x${2}x${3}",
		      }
		      );
	die "ERROR: Need at least 1 dimension replication > 1!\n" if ($1 == $2 and $2 == $3 and $3 == 1);
    }

    if (! defined($saveFile)) {
	$saveFile = basename ($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_" . $replicate->{STRING} . ".bgf";
    }
    $doCenter = 0 if (! defined($doCenter) or $doCenter !~ /1|yes/i);
    $bondImages = 1 if (! defined($bondImages) or $bondImages !~ /0|no/i);
    print "Done\n";
}
