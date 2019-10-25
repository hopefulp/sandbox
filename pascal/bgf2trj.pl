#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use Packages::General qw(FileTester);
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::AMBER qw(CreateAmberTrj);
use Packages::BOX qw(GetBox);
use Packages::LAMMPS qw(CreateLAMMPSTrj);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub createTrj;
sub init;
sub createLAMMPSHeaders;

my ($FILES, $save, $saveFunc, $saveType);

$|++;

&init;
open OUTFILE, "> $save" or die "ERROR: Cannot create $saveType Trj file $save: $!\n";
print OUTFILE "TITLE: $saveType Trajectory file created by bgf2trj.pl on " . scalar(localtime(time())) . "\n";
&createTrj($FILES, $saveFunc, \*OUTFILE);
print "Created $saveType trajectory $save\n";
close OUTFILE;

sub createTrj {
    my ($flist, $trjWriter, $outHandle) = @_;
    my ($file, $tot, $ATOMS, $BONDS, $HEADERS, $count, $BOX);

    $count = 0;
    for $file (@{ $flist }) {
	print "Reading BGF file $file...\r";
	($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($file, 1,0);
	$BOX = GetBox($ATOMS, undef, $HEADERS);
	$BOX = createLammpsHeaders($BOX, $count, scalar(keys %{ $ATOMS}));
	$trjWriter->($ATOMS, $BOX, $outHandle);
	$count++;
    }
    print "Reading BGF file...Read " . scalar(@{ $flist }) . " file...Done                                          \n";
}

sub createLammpsHeaders {
    my ($box, $tStep, $numAtms) = @_;
    my (%DATA);

    $DATA{"TIMESTEP"}[0] = $tStep;
    $DATA{"NUMBER OF ATOMS"}[0] = $numAtms;
    $DATA{"BOX BOUNDS"}[0]{lo} = $box->{X}{lo};
    $DATA{"BOX BOUNDS"}[0]{hi} = $box->{X}{hi};
    $DATA{"BOX BOUNDS"}[1]{lo} = $box->{Y}{lo};
    $DATA{"BOX BOUNDS"}[1]{hi} = $box->{Y}{hi};
    $DATA{"BOX BOUNDS"}[2]{lo} = $box->{Z}{lo};
    $DATA{"BOX BOUNDS"}[2]{hi} = $box->{Z}{hi};

    return \%DATA;
}

sub init {
    my (%OPTS);

    getopt('bst',\%OPTS);
    die "usage: $0 -b \"bgf file(s) or location\" -t (type=lammps(default) or amber) -s (savename)\n"
        if (! exists($OPTS{b}));
    print "Initializing...";
    open FINDCMD, "ls $OPTS{b} | grep \"\.bgf\$\" |" or die "ERROR: Cannot locate files in \"$OPTS{b}\"\n";
    while (<FINDCMD>) {
        chomp;
        push @{ $FILES }, $_;
    }
    close FINDCMD;
    die "ERROR: No valid .bgf file found while searching \"$OPTS{b}\"!\n"
        if (! $FILES);
    if (! exists($OPTS{s})) {
        $save = basename($FILES->[0]);
        $save =~ s/\.\w+$//;
        $save .= "_orientation.dat";
    } else {
	$save = $OPTS{s};
    }
    $saveType = "lammps" if (! exists($OPTS{t}) or $OPTS{t} !~ /lammps/i);
    $saveType = "amber" if (exists($OPTS{t}) and $OPTS{t} !~ /lammps/i);
    $saveFunc = \&CreateLAMMPSTrj;
    $saveFunc = \&CreateAmberTrj if ($saveType ne "lammps");
    $saveType = uc $saveType;
    print "Done\n";
}
