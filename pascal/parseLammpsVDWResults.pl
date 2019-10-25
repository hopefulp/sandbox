#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub getLammpsVals;
sub writeCSVfile;
sub numerically { ($a<=>$b); }

my ($csvFile, $engFiles);
my ($DATA);

$|++;
&init;
print "Getting energy from jag files in '$engFiles'...";
$DATA = getLammpsVals($engFiles);
print "Done\n";
&writeCSVfile($DATA, $csvFile);
print "Done\n";

sub writeCSVfile {
    my ($engData, $outFile) = @_;
    my ($i, $j, $k, @headers, @dist, $single);

    if (! defined($outFile)) {
	for $i (keys %{ $engData }) {
	    $outFile .= "${i}_";
	}
	chop $outFile;
	$outFile .= "_results.csv";
    }
    print "Writing data $outFile...";

    open OUTDATA, "> $outFile" or die "ERROR: Cannot write to $outFile: $!\n";
    for $i (keys %{ $engData }) {
	print OUTDATA "${i},\n";
	@dist = sort numerically  keys %{ $engData->{$i} };
	@headers = keys %{ $engData->{$i}{ $dist[0] } };
	print OUTDATA "dist,";
	for $j (@headers) {
	    print OUTDATA "${j},";
	}
	print OUTDATA "\n";
	for $j (@dist) {
	    print OUTDATA "${j},";
	    for $k (@headers) {
		$engData->{$i}{$j}{$k} = "" if (! exists($engData->{$i}{$j}{$k}));
		print OUTDATA "$engData->{$i}{$j}{$k},";
	    }
	    print OUTDATA "\n";
	}
	print OUTDATA "\n";
    }
    close OUTDATA;
}

sub getLammpsVals {
    my ($fileLoc) = $_[0];
    my ($cmd, $ENG, $name, $dist, $type, $val, $files, $i);
    
    $cmd = "ls $fileLoc";
    open FILELIST, "$cmd |" or die "ERROR: No valid files found while searching '$cmd':$!\n";
    while (<FILELIST>) {
	chomp;
	push @{ $files }, $_;
    }
    close FILELIST;

    for $i (@{ $files }) {
	$name = basename($i);
        if ($name =~ /^(\D+)_(\d+\.?\d*)_(\S+)\.?\w*$/) {
	    ($name, $dist, $type) = ($1, $2, $3);
            $type =~ s/graphite_//;
            $type =~ s/_log\.?\w*$//;
            $name =~ s/_mul//g;
 	    $name =~ s/_esp//g;
	    undef($val);
	    open ENGFILE, $i or die "ERROR: Cannot open $i: $!\n";
            while (<ENGFILE>) {
		chomp;
		if ($_ =~ /E_vdwl\s+\=\s+(\-?\d+\.\d+)/) {
		    $val = $1;
		}
            }
            close ENGFILE;
	    $ENG->{$name}{$dist}{$type} = $val if (defined($val));
        }
    }
    die "ERROR: Nothing matched!\n" if (! defined($ENG));
    return $ENG;
}

sub init {
    my (%OPTS);

    getopt('ls',\%OPTS);
    die "usage:$0 -l lammps log file(s) -s (save file)\n" if (! exists($OPTS{l}));
    print "Initializing...";
    die "ERROR: No valid files found while searching '$OPTS{l}'\n" if(system("ls $OPTS{l} > /dev/null"));
    ($engFiles, $csvFile) = ($OPTS{l}, $OPTS{s});
    print "Done\n";
}
