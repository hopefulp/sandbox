#!/usr/bin/perl -w
BEGIN {
	unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::General qw(FileTester);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub createLeapFile;
sub runLeap;

my ($leaprc, $pdbFile, $prefix);

$|++;
&init;
print "Creating leaprc file using $leaprc...";
&createLeapFile($leaprc, $prefix, $pdbFile);
print "Done\nCreating AMBER $prefix files from $pdbFile...";
&runLeap($prefix);
system("rm -fr leaprc");
print "Done\n";

sub runLeap {
    my ($mol) = $_[0];
    if (system "/exec/amber9/exe/tleap >& ${mol}.out") {
	print "Error Occurred. Email ${mol}.out file to tpascal\@wag.caltech.edu. Here is what it contains:\n";
	system "cat ${mol}.out";
    }
}

sub createLeapFile {
    my ($leapFile, $molname, $pdbfile) = @_;

    die "ERROR: Cannot copy $leapFile to current directory!\n" if (system("cp $leapFile ./leaprc"));
    open MYLEAPRC, ">> leaprc" or die "Cannot write to leaprc: $!\n";
    print MYLEAPRC "gaff = loadamberparams gaff.dat\n";
    print MYLEAPRC "prot = loadpdb $pdbfile\n";
    print MYLEAPRC "saveamberparm prot $molname" . ".prmtop $molname" . ".inpcrd\n";
    print MYLEAPRC "quit\n";
    close MYLEAPRC;
}
    
sub init {
    my (%OPTS);

    getopt('lps',\%OPTS);
    die "usage: $0 -p pdb file -s (save prefix) -l (leaprc location)\n" if (! exists($OPTS{p}));
    print "Initializing...";
    ($pdbFile, $leaprc, $prefix) = ($OPTS{p}, $OPTS{l}, $OPTS{s});
    FileTester($pdbFile);
    $leaprc = "/home/yjn1818/amber/leaprc" if (! defined($leaprc));
    if (! -e $leaprc or ! -r $leaprc or ! -T $leaprc) {
	print "invalid amber leaprc file... using default...";
        $leaprc = "/home/yjn1818/amber/leaprc";
    }
    $prefix = basename($pdbFile) if (! defined($prefix));
    $prefix =~ s/\.\w+$//;
    print "Done\n";
}

