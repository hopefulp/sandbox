#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::MolData;
use Getopt::Std qw(getopt);

sub init;
sub getCharges;
sub updateCharges;

my ($molStructure, $respChgFile, $saveName, $fileType, $fileName);
my ($CHARGES);

$|++;
&init;

print "Gettting data from $fileType file $fileName...";
$molStructure->read($fileName, $fileType);
print "Done\nParsing charge file $respChgFile....";
$CHARGES = getCharges($respChgFile);
print "Done\nUpdating charges...";
&updateCharges($molStructure, $CHARGES);
print "Done\nWriting to $fileType file $saveName...";
$molStructure->write($saveName, $fileType);
print "Done\n";

sub updateCharges {
    my ($molData, $chargeData) = @_;
    my ($tot, $i, $atomTot);

    $tot = scalar(keys %{ $chargeData });
    $atomTot = $molData->count->{"atoms"};
    die "ERROR: Unequal number of atoms in data file ($atomTot) vs charges ($tot)!\n"
	if ($tot != $atomTot);

    for $i (1 .. $tot) {
	$molData->atoms($i)->("charge", $chargeData->{$i});
    }
}



sub getCharges {
    my ($infile) = $_[0];
    my (%DATA, $count);
    
    $count = 0;
    open CHRG, $infile or die "ERROR: Cannot open resp charge file $infile: $!\n";
    while (<CHRG>) {
	chomp;
	if ($_ =~ /^\s*\-?\d+\.\d+/) {
	    while ($_ =~ /(\-?\d+\.\d+)/g) {
		$count++;
		$DATA{$count} = $1;
	    }
	}
    }
    close CHRG;

    die "ERROR: No valid data while reading $infile!\n" if (! $count);

    return \%DATA;
}

sub init {
    my (%OPTS);

    getopt('ftsc', \%OPTS);
    for ("f", "c") {
	die "usage: $0 -f file name -c charge file -t (file type) -s (savename)\n" 
	    if (! exists($OPTS{$_}));
    }

    print "Initialzing...";
    $molStructure =  Packages::MolData->new();
   ($fileName, $fileType, $saveName, $respChgFile) = ($OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{c});
    $molStructure->testFile($fileName);
    $molStructure->testFile($respChgFile);
    if (! defined($fileType)) {
        $fileType = "bgf";
        if ($fileName =~ /\.(\w+)$/) {
            $fileType = lc $1;
        }
    }

    $saveName = $molStructure->getFileName($fileName) if (! defined($saveName));
    print "Done\n";
}
