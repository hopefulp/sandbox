#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts/";
}

use strict;
use Packages::MolData;
use Getopt::Std qw(getopt);

sub init;

my ($molStructure, $fileName, $fileType, $readFunc);

$|++;
&init;
print "Gettting data from $fileType file $fileName...";
$molStructure->read($fileName, $fileType);

print "Done\n";

print "";


sub init {
    my (%OPTS);

    getopt('ft', \%OPTS);
    die "usage: $0 -f file name -t (file type)\n" if (! exists($OPTS{f}));

    print "Initialzing...";
    $molStructure =  Packages::MolData->new();
   ($fileName, $fileType) = ($OPTS{f}, $OPTS{t});
    $molStructure->testFile($fileName);
    if (! defined($fileType)) {
	$fileType = "bgf";
	if ($fileName =~ /\.(\w+)$/) {
	    $fileType = lc $1;
	}
    }

    print "Done\n";
}
