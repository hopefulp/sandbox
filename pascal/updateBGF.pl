#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF);
use Packages::General qw(FileTester Trim);
use Getopt::Std qw(getopt);

sub main;
sub init;
sub updateBGF;

my ($refBGF, $modBGF, $saveName, $field, $newBox);

$|++;
&main;

sub main {
    my ($refATOMS, $refBONDS, $HEADERS, $modATOMS, $tmp, $FIELDS, $HEADERS1);

    $FIELDS = &init;
    print "Obtaining reference BGF Info from $refBGF...";
    ($refATOMS, $refBONDS, $HEADERS) = GetBGFFileInfo($refBGF, 1);
    print "Done\nObtaining coordiantes/charges from $modBGF...";
    ($modATOMS, $tmp, $HEADERS1) = GetBGFFileInfo($modBGF, 1);
    print "Done\nUpdating \"$field\"...";
    updateBGF($refATOMS, $modATOMS, $FIELDS);
    print "Done\nCreating BGF file $saveName...";
    if ($newBox) {
	addHeader($refATOMS, $HEADERS1);
    } else {
	addHeader($refATOMS, $HEADERS);
    }
    createBGF($refATOMS, $refBONDS, $saveName);
    print "Done\n";
}

sub updateBGF {
    my ($ref, $mod, $fieldList) = @_;
    my ($atomC, $modName, $atom, $field, $refNames);

    for $atomC (keys %{ $ref }) {
	$refNames->{ Trim($ref->{$atomC}{"ATMNAME"}) } = $atomC;
    }

    for $atomC (keys %{ $mod }) {
	die "ERROR: Atom # $atomC ($mod->{$atomC}{ATMNAME}) is not present in reference file\n"
	    if (! exists($ref->{$atomC}));
	$modName = Trim($mod->{$atomC}{"ATMNAME"});
	die "ERROR: Atom name \"$modName\" not found in reference file!\n" if (! exists($refNames->{$modName}));
	$atom = $ref->{ $refNames->{$modName} };
	for $field (@{ $fieldList }) {
	    next if (! exists($atom->{$field}));
	    $atom->{$field} = $mod->{$atomC}{$field};
	}
    }
}
	
sub init {
    my (%OPTS, @FIELDS, $i);

    getopt('brts',\%OPTS);
    ($refBGF, $modBGF, $field, $saveName) = ($OPTS{b},$OPTS{r},$OPTS{t},$OPTS{s});
    for ($refBGF, $modBGF, $field) {
	die "usage: $0 -b updateBGF -r referenceBGF -t \"fields\" -s [saveBGF]\n"
	    if (! defined($_));
    }
    
    print "Initializing...";
    FileTester($refBGF);
    FileTester($modBGF);

    if (! $saveName) {
	$saveName = $refBGF;
	$saveName =~ s/\.\w+$/_mod\.bgf/;
    }

    if ($field =~ /\s+/) {
	@FIELDS = split /\s+/, $field;
    } else {
	$FIELDS[0] = $field;
    }

    $newBox = 0;
    for $i (@FIELDS) {
	if (lc($i) eq "box") {
	     $newBox =1;
	     last;
	}
    }
    print "Done\n";
    return \@FIELDS;
}
    
