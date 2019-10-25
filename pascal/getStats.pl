#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;

use Getopt::Std qw(getopt);
use File::Basename qw(basename);

use Packages::General qw(FileTester STDev GetStats GetEquilPoint);

sub init;
sub parseDataFile;
sub getBestVal;
sub numerically { ($a<=>$b); }

my ($DATA, $dataFile, $STATS, $tolerance, $findex);
my ($equilDATA, $isFound);

$|++;
&init;
print "Getting data from $dataFile...";
$DATA = parseDataFile($dataFile, $findex);
print "Done\nDetermining Equilibration point...";
($isFound, $equilDATA) = GetEquilPoint($DATA->{VALS}, $tolerance);
print "not found...will use all data..." if (! $isFound);
print "Done\nCollecting Stats...";
$STATS = GetStats($equilDATA);
print "Done\n";
print "EQUIL PT: $isFound\n" if ($isFound);
print "AVG: $STATS->{AVG}\nSTDEV: $STATS->{STDEV}\n";
&getBestVal($STATS,$equilDATA, $DATA) if ($isFound);

sub getBestVal {
    my ($stats, $data, $origdata) = @_;
    my ($residuals, $i, $bestval, @tmp);

    for $i (0 .. $#{ $data }) {
	$residuals->{ abs($stats->{AVG}-$data->[$i]) } = $i;
    }

    @tmp = sort numerically keys %{ $residuals };
    $i = shift @tmp;
    $bestval = $data->[ $residuals->{ $i } ];
    print "BEST: $bestval [diff: $i";
    if(defined($origdata)) {
	print " #:";
	for $i (keys %{ $origdata->{VALS} }) {
	    if($origdata->{VALS}{$i} == $bestval) {
		print " $i ( $origdata->{COUNT}{$i} )";
	    }
	}
    }
    print " ]\n";
}

sub parseDataFile {
    my ($inFile, $field) = @_;
    my (%data, $counter, @tmp);

    $counter = 0;
    open DATFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
    while (<DATFILE>) {
	chomp;
	last if ($_ =~ /^##break/);
	if ($_ =~ /^\s*(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*.*$)/) {
	    @tmp = split /\s+/,$2;
	    if($field > $#tmp) {
		close DATFILE;
		die "ERROR: Field number " . ($field+1) . " does not exists!\n";
	    }
	    $data{VALS}{$1} = $tmp[$findex];
	    $data{COUNT}{$1} = $counter+1;
	    $counter++;
	} elsif ($_ =~ /^\s*(\-?\d+\.?\d*)\s*$/) {
	    $data{VALS}{$counter} = $1;
	    $data{COUNT}{$counter} = $counter+1;
	    $counter++;
	}
    }
    close DATFILE;
    
    die "ERROR: $inFile does not contain any valid information!\n" if (! %data);
    return \%data;
}

sub init {
    my (%OPTS);
    getopt('dti',\%OPTS);

    die "usage: $0 -d data file -t (tolerance %) -i (field index=2)\n" if (! exists($OPTS{d}));
    print "Initializing...";
    ($dataFile, $tolerance, $findex) = ($OPTS{d}, $OPTS{t}, $OPTS{i});
    FileTester($dataFile);
    $tolerance = 1 if (! defined($tolerance));
    die "ERROR: Expected integer/decimal for tolerance. Got \"$tolerance\"\n"
	if ($tolerance !~ /^\d+\.?\d*/);
    print "warning.. tolerance > 100% of average!..." if ($tolerance > 100);
    $tolerance /= 100;
    $findex = 2 if (! defined($findex) or $findex !~ /^\d+$/ or $findex<1);
    $findex -= 2;
    print "Done\n";
}
