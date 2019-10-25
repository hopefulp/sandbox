#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use warnings;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use Packages::General qw(FileTester);
use Packages::ManipAtoms qw(SplitAtomsByMol SelectAtoms BuildAtomSelectionString GetMols GroupAtomsByField);

sub init;
sub showUsage;

my ($bgfFile, $saveFile, $SELECT, $selection, $field);
my ($ATOMS, $BONDS, $HEADERS, $MOLS);

$|++;
&init;
print "Parsing bgf file...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
if(defined($selection)) {
  print "selecting atoms...";
  $SELECT = SelectAtoms($selection, $ATOMS);
  die "No selected atoms found in bgf file\n" if (! keys %{ $SELECT });
}
print "Done\nGrouping atoms based on $field...";
&GroupAtomsByField($ATOMS, $BONDS, $field, $SELECT);
print "Done\nCreating bgf file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub init {
    my (%OPTS, $atomSelect, $usage);

    &getopt('bafs',\%OPTS);
    $usage = showUsage();
    die "$usage\n" if (! exists($OPTS{b}));

    print "Initializing...";
    ($bgfFile, $saveFile, $atomSelect, $field) = ($OPTS{b}, $OPTS{s}, $OPTS{a}, $OPTS{f});
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+//;
	$saveFile .= "_mod.bgf";
    }
    
    $field = "CHAIN" if (! defined($field) or $field !~ /^CHAIN|RESNUM|RESID|RES|MOL|MOLECULE|MOLID|MOLECULEID|MOLSIZE|NUMBONDS|FFTYPE$/i);
    $field = uc($field);
    $field = "RESNUM" if ($field =~ /RESID/);
    $field = "MOLECULEID" if ($field =~ /MOL/ && $field !~ /SIZE/);
    $field = "MOLSIZE" if ($field =~ /SIZE/);
    if (defined($atomSelect)) {
	$selection = BuildAtomSelectionString($atomSelect);
    }
    print "Done\n";
}

sub showUsage {
    my ($usage) = <<DATA;
usage: $0 -b bgf_file -a (atom_selection) -f (field) -s (save_name)
Options:
	-b bgf_file (Required): location of bgf file
	-a atom_selection (Optional): any valid field expression. E.g. "resname eq 'WAT' and zcoord > 50"
	-f field (Optional): chain|resid|resname|molecule(id)|molsize|numbonds. Default: chain
	-s save_name (Optional): name of new bgf file
DATA

    return $usage;
}
