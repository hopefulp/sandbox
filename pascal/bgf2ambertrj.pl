#!/usr/bin/perl
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::General qw(FileTester);
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::AMBER qw(CreateAmberTrj);
use Packages::BOX qw(GetBox);

sub main;
sub init;
sub numerically;

$|++;

die "usage: $0 bgf_file|directory amber_traj_name\n"
    if (! @ARGV or $#ARGV < 1);
my ($bgfLoc, $saveName) = @ARGV;
my ($scaled) = 0;

&main;

sub main {
    my ($FLIST, $file);
    my ($ATOMS, $BONDS, $BOX, $HEADERS);

    $FLIST = &init;
    print "Creating AMBER trajectory file $saveName\n";
    open OUTFILE, "> $saveName" or die "ERROR: Cannot create Amber Trj file $saveName: $!\n";
    print OUTFILE "TITLE: AMBER Trajectory file created by bgf2ambertrj.pl on " .
	scalar(localtime(time())) . "\n";
    for $file (@{ $FLIST }) {
	print "Reading BGF file $file...";
	($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($file, 1);
	$BOX = GetBox($ATOMS, undef, $HEADERS);
	CreateAmberTrj($ATOMS, $BOX, \*OUTFILE);
        print "Done\r";
    }
    close OUTFILE or die "ERROR: Cannot finalize file $saveName: $!\n";
    print "\nDone\n";
}

sub init {
    my (@FILES, $findCmd);

    $findCmd = "find " . dirname($bgfLoc) . " -name '*.bgf' -print" if (! -e $bgfLoc);
    $findCmd = "ls $bgfLoc" if (-e $bgfLoc);
    die "ERROR: No valid file found while searching $bgfLoc!\n"
        if (! open(DATFILES, "$findCmd |"));
    while (<DATFILES>) {
        chomp;
        if (-e $_ and -r $_ and -T $_) {
            push @FILES, $_;
        }
    }
    close DATFILES;
    die "ERROR: No valid files found in path \"$bgfLoc\"\n"
        if (! @FILES);
    print "Found " . ($#FILES + 1) . " BGF files\n";

    return \@FILES;
}

sub numerically {
    ($a<=>$b);
}
