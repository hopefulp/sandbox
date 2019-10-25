#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::AMBER qw(getTopInfo ParseAmberTrj getOpts getConn GetAmberByteOffset);
use Packages::General qw(FileTester TrjSelections);
use Packages::FileFormats qw(createBGF createHeaders addHeader);
use Packages::BOX qw(MakeBox);
use File::Basename;

sub init;
sub numerically;
sub createFile;

die "usage: $0 topFile trjFile \"trajectory selection\" [savePrefix]\n"
    if (! @ARGV || $#ARGV < 3);

my ($topFile, $trjFile, $selection, $savePrefix) = @ARGV;
my ($OPTS, $SELECT, $DATA, $totAtms);
my ($printStr, $VEC, $link, $atmCounter, $trjName);

$|++;
print "Initializing...";
&init;
print "Done\nParsing AMBER topology file $topFile...";
($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
print "Done\nGenerating Connectivities...";
$DATA->{BONDS} = getConn(\%{ $DATA->{"BONDLIST"} }, $DATA->{"ATOMS"});
print "Done\n";
$printStr = "Parsing AMBER trajectory $trjFile...";
&GetAmberByteOffset($SELECT, $trjFile, $totAtms);
ParseAmberTrj($DATA->{ATOMS}, $trjFile, $SELECT, $totAtms, \&createFile, $printStr, $DATA->{BONDS});

sub init {
    FileTester($topFile);
    FileTester($trjFile);

    $OPTS = &getOpts;
    $savePrefix = basename($trjFile) if (! defined($savePrefix));
    $savePrefix =~ s/\.\w+$//;
 
    if (-e $selection) {
	open SEL, $selection or die "ERROR: Cannot open file $selection: $!\n";
	$selection = "";
	while (<SEL>) {
	    chomp;
	    $selection .= "$_ ";
	}
	close SEL;
    }   
    $SELECT = TrjSelections($selection);
}

sub createFile {
    my ($ATOMS, $BBOX, $frameNum, $BONDS) = @_;
    my ($i, $BOX, $counter, @tmp);
    my ($bgfFile) = "${savePrefix}_${frameNum}.bgf";

    $BOX->{1}{DATA} = 90;
    $counter = 2;
    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (0 .. $#tmp) {
	$BOX->{$counter}{DATA} = $BBOX->{($i + 2)}{DATA};
	$counter++;
    }
    $BOX = MakeBox($BOX);
    my ($HEADERS) = createHeaders($BOX, $bgfFile);
    addHeader($ATOMS, $HEADERS);
    createBGF($ATOMS, $BONDS, $bgfFile);
}
