#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::General qw(FileTester TrjSelections);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF createHeaders);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset ConvertLammpsBox);
use Packages::ManipAtoms qw(GetMols);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub storeIndices;
sub init;
sub createGrpFile;
sub getTrjData;
sub writeBGF;
sub numerically { ($a<=>$b); }

my ($trjs, $bgfFile, $saveName);
my ($ATOMS, $BONDS, $DATA, $HEADERS);

$|++;
&init;

print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 1);
&GetMols($ATOMS, $BONDS);
print "Done\n";
&getTrjData($trjs, $ATOMS, \%{ $DATA });
print "Creating Group file ${saveName}...";
&createGrpFile($DATA, $saveName);
print "Done\nCreating shell bgf...";
&writeBGF($ATOMS, $BONDS, $HEADERS, $saveName);
print "Done\n";

sub writeBGF {
    my ($atoms, $bonds, $headers, $save) = @_;
    my ($i);

    $save =~ s/\.\w+//g;
    $save .= "_shells.bgf";

    for $i (keys %{ $atoms }) {
	$atoms->{$i}{CHAIN} = uc(chr(65 + $atoms->{$i}{ASSIGNED}));
    }
    &addHeader($atoms, $headers);
    &createBGF($atoms, $bonds, $save);
}

sub createGrpFile {
    my ($data, $saveFile) = @_;
    my ($i, @tmp, @tmp2, $j);
    my ($groups, $start, $prev, $counter);

    $i = 0;
    @tmp = sort numerically keys %{ $data };
    while ($i < $#tmp) {
	@tmp2 = keys %{ $data->{$tmp[$i]} };
	if (! @tmp2) {
	    splice @tmp, $i, 1;
	    next;
	}
	$i++;
    }
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
                $groups .= "${start} - ${prev} ";
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

    close OUTDATA;
}

sub storeIndices {
    my ($indices, $atoms, $groupID) = @_;
    my ($i, $molList, $j, $tStep, $BOX);

    $BOX = ConvertLammpsBox($indices->{"BOX BOUNDS"});
    $HEADERS = createHeaders($BOX,  $indices->{TIMESTEP}[0]);
    for $i (keys %{ $indices->{ATOMS} }) {
	next if exists($atoms->{$i}{ASSIGNED});
	$molList = $atoms->{$i}{MOLECULE}{MEMBERS};
	for $j (keys %{ $molList }) {
	    $atoms->{$j}{ASSIGNED} = $groupID;
	    $DATA->{$groupID}{$j} = 1;
	}
    }
}

sub getTrjData {
    my ($list, $atoms, $data) = @_;
    my ($i, $pStr, $SELECT, $LAMMPSOPTS);

    $LAMMPSOPTS->{imaged} = 1;
    $LAMMPSOPTS->{scaled} = 1;

    $SELECT = TrjSelections("1");
    for $i (0 .. $#{$list}) {
	&GetLammpsByteOffset($SELECT, $list->[$i], undef);
	&ParseLAMMPSTrj($atoms, $list->[$i], $SELECT, "atom", \&storeIndices, undef, $i);
    }
}

sub init {
    my (%OPTS, $trjStr);

    getopt('lbs',\%OPTS);

    ($trjStr, $bgfFile,$saveName) = ($OPTS{l},$OPTS{b},$OPTS{s});

    for ($trjStr, $bgfFile) {
	die "usage:$0 -b bgf file -l lammps trj(s) -s (savename)\n" if (! defined($_));
    }
    print "Initializing...";

    FileTester($bgfFile);
    while ($trjStr =~ /(\S+)/g) {
	push @{ $trjs }, $1 if (-e $1);
    }
    die "ERROR: No valid Lammps trajectories found while searching \"$trjStr\"!\n"
	if (! defined($trjs));
    if (! defined($saveName)) {
	$saveName = basename($trjs->[0]);
	$saveName =~ s/\.\w+$//;
	$saveName .= ".grps";
    }
    print "Done\n";
}
