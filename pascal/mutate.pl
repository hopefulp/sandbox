#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts/";
}

use strict;
use Packages::General qw(FileTester GetSelections CoM CrossProduct);
use Packages::FileFormats qw(GetBGFFileInfo addHeader createBGF sortByRes GetBGFAtoms);
use Packages::ManipAtoms qw(GetAtmList);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::Math::MatrixReal;

sub init;
sub mutateResidues;
sub GetResList;
sub numerically { ($a<=>$b); }
sub matchAtoms;
sub InsertResidue;
sub growAtom;
sub buildLinkedList;
sub buildHashTable;
sub deleteAtom;
sub findNode;

my ($bgfFile, $resFile, $saveFile, $resList);
my ($ATOMS, $BONDS, $HEADERS, $resATOMS, $resBONDS);
my ($selection, $RES, $aList, $mList);
$|++;
&init;
print "Parsing bgf file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nParsing residue file $resFile...";
($resATOMS, $resBONDS) = GetBGFFileInfo($resFile, 0);
print "Done\nGetting mutation residue information...";
$RES = GetSelections($selection, 0);
$RES = GetAtmList($RES, $ATOMS);
$RES = GetResList($ATOMS, $RES);
$aList = buildLinkedList($ATOMS, $BONDS);
$mList = buildLinkedList($resATOMS, $resBONDS);
die "ERROR: No valid atoms selected!\n" if (! keys %{ $RES });
print "Done\n";
&mutateResidues($ATOMS, $aList, $RES, $resATOMS, $resBONDS, $mList);
print "Creating bgf file $saveFile...";
($ATOMS, $BONDS) = buildHashTable($aList);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub buildHashTable {
    my ($linkedList) = @_;
    my (%atoms, %bonds, $node, $head, $index, $bond, @bondList, %bondCheck, @tmp);

    $node = $linkedList->{HEAD};
    $head = $node;
    $node->{data}{INDEX} = 1;
    $node->{data}{FLAG} = 1;
    $node = $head->{next};
    $index = 2;
    while (! exists($node->{data}{FLAG})) {
	$node->{data}{INDEX} = $index;
	$node->{data}{FLAG} = 1;
	$index++;
	$node = $node->{next};
    }

    $node = $head;
    while ($node->{data}{FLAG} == 1) {
	$index = $node->{data}{INDEX};
	@tmp = grep !/^BOND|MOL|NODE/, keys %{ $node->{data} };
	for (@tmp) {
	    $atoms{$index}{$_} = $node->{data}{$_};
	}
	delete $atoms{$index}{NODE};
	@bondList = grep /^BOND/, keys %{ $node };
	%bondCheck = ();
	for $bond (@bondList) {
	    if (keys %{ $node->{$bond} } and keys %{ $node->{$bond}{data} } and 
		! exists($bondCheck{$node->{$bond}{data}{INDEX}})) {
		push @{ $bonds{$index} }, $node->{$bond}{data}{INDEX};
		$bondCheck{$node->{$bond}{data}{INDEX}} = 1;
	    }
	}
	$node->{data}{FLAG} = 2;
	$node = $node->{next};
    }

    return (\%atoms, \%bonds);
}

sub buildLinkedList {
    my ($atoms, $bonds) = @_;
    my ($root, $i, @indices, $last, $node, $prev, %atomList, $j, @bondList, $k);

    @indices = sort numerically keys %{ $atoms };
    $i = shift @indices;
    $root = {data => $atoms->{$i}, next => undef, prev => undef};
    $last = $root;
    $prev = undef;
    $atoms->{$i}{NODE} = $root;
    @bondList = grep /^BOND/, keys %{ $atoms->{$i} };
    for $j (@bondList) {
	$root->{$j} = $atoms->{$i}{$j};
	delete $atoms->{$i}{$j};
    }

    for $i (@indices) {
	$node = {data => $atoms->{$i}, next => undef, prev => $last};
	$last->{next} = $node;
	$prev = $last;
	$last = $node;
	$atoms->{$i}{NODE} = $node;
	@bondList = grep /^BOND/, keys %{ $atoms->{$i} };
	for $j (@bondList) {
	    $node->{$j} = $atoms->{$i}{$j};
	    delete $atoms->{$i}{$j};
	}
    }
    $root->{prev} = $last;
    $last->{next} = $root;
    $root->{ISHEAD} = 1;
    $last->{ISTAIL} = 1;

    for $i (keys %{ $bonds }) {
	$node = $atoms->{$i}{NODE};
	for $j (0 .. $#{ $bonds->{$i} }) {
	    $k = $j;
	    $k++ while (exists($node->{"BOND" . $k}));
	    $node->{"BOND" . $k} = $atoms->{ $bonds->{$i}[$j] }{NODE};
	}
    }

    $atomList{HEAD} = $root;
    $atomList{TAIL} = $last;
    $atomList{COUNT} = scalar(@indices) + 1;
    return \%atomList;
}

sub GetResList {
    my ($atoms, $res) = @_;
    my (%RES, $i);

    for $i (keys %{ $res }) {
	die "ERROR: Atom $i not found!\n" if (! exists($atoms->{$i}));
	$RES{ $atoms->{$i}{RESNUM} }{$i} = 1;
    }
    return \%RES;
}

sub mutateResidues {
    my ($atoms, $atomsList, $resInfo, $mutateAtoms, $mutateBonds, $mutateList) = @_;
    my ($i, $resNum, @tmp, $mutateName, $printStr, $count, $j, $resName);
    my ($mutateCount, $offset, $startAtm, $resAtoms, $resList, $mutateStr);
    
    $mutateCount = $mutateList->{COUNT};
    $mutateName = $mutateList->{HEAD}{data}{RESNAME};

    for $i (sort numerically keys %{ $resInfo }) {
	@tmp = sort numerically keys %{ $resInfo->{$i} };
	$offset = (scalar(@tmp) - $mutateCount);
	$startAtm = $tmp[0];
	$resNum = $atoms->{ $tmp[0] }{RESNUM};
	$resName = $atoms->{ $tmp[0] }{RESNAME};
	$printStr = "Mutating residue # $i ($resName) to $mutateName...";
	    print "${printStr}\r";
	$mutateStr .= "$resName (#$resNum) ";
	$resAtoms = matchAtoms($atomsList, $resInfo->{$i}, $mutateList, $resNum);
	next if (! $resAtoms);
	$resList = buildLinkedList($resAtoms, undef);
	for $j (@tmp) {
	    $startAtm = &deleteAtom($atoms->{$j}{NODE}, \%{ $atomsList });
	}
	if (keys %{ $startAtm }) {
	    $offset = $startAtm->{next};
	    $startAtm->{next} = $resList->{HEAD};
	    $resList->{HEAD}{prev} = $startAtm;
	    $resList->{TAIL}{next} = $offset;
	    $offset->{prev} = $resList->{TAIL};
	    $atomsList->{COUNT} += $resList->{COUNT};
	    if (exists($offset->{ISHEAD}) and ! exists($startAtm->{ISTAIL})) {
		$atomsList->{HEAD} = $resList->{HEAD};
		delete $offset->{ISHEAD};
	    } else {
		delete $resList->{HEAD}{ISHEAD};
	    }
	    if (exists($startAtm->{ISTAIL})) {
		$atomsList->{TAIL} = $resList->{TAIL};
		delete $startAtm->{ISTAIL};
	    } else {
		delete $resList->{TAIL}{ISTAIL};
	    }
	} else {
	    for (keys %{ $resList }) {
		delete $atomsList->{$_};
		$atomsList->{$_} = $resList->{$_};
	    }
	}
	$count++;
    }
    die "Mutating residues: ERROR: No residues mutated!\n" if (! $count);
    printf "%-" . (length($printStr) + 10) . "s\n", "Mutating residues: Sucessfully mutated $mutateStr -> $mutateName";
}

sub printList {
    my ($aList) = @_;
    my ($node, $i, @tmp);
    
    return () if ($aList->{COUNT} < 1);
    $node = $aList->{HEAD};
    print "... HEAD $node->{data}{ATMNAME}: (# $node->{data}{INDEX})\n";
    @tmp = grep /^BOND/, keys %{ $node };
    for $i (@tmp) {
	print "\t\t";
	next if (! keys %{ $node->{$i}{data} });
	print "$i: $node->{$i}{data}{ATMNAME} ($node->{$i}{data}{INDEX}) ";
    }
    print "\n";
    $node = $node->{next};
    while (! exists($node->{ISHEAD})) {
	print "\t $node->{data}{ATMNAME}: (# $node->{data}{INDEX})\n";
	@tmp = grep /^BOND/, keys %{ $node };
	print "\t\t";
	for $i (@tmp) {
	    next if (! keys %{ $node->{$i}{data} });
	    print "$i: $node->{$i}{data}{ATMNAME} ($node->{$i}{data}{INDEX}) ";
	}
	print "\n";
	$node = $node->{next};
    }
}
    
sub deleteAtom {
    my ($node, $atomsList) = @_;
    my ($prev, $next, $head, $i, $returnNode);
    
    $prev = $node->{prev};
    $next = $node->{next};
    $returnNode = $prev;

    if (exists($node->{ISHEAD})) {
	$atomsList->{HEAD} = $next;
	$next->{ISHEAD} = 1;
    } elsif (exists($node->{ISTAIL})) {
	$atomsList->{TAIL} = $prev;
	$prev->{ISTAIL} = 1;
	#$returnNode = $prev->{prev};
    }
    $prev->{next} = $next;
    $next->{prev} = $prev;
    for $i (keys %{ $node->{data} }) {
	delete $node->{data}{$i};
    }
    for $i (keys %{ $node }) {
	delete $node->{$i};
    }
    undef($node);
    $atomsList->{COUNT}--;
    return $returnNode;
}

sub matchAtoms {
    my ($atomsList, $resAtomList, $resList, $resNum) = @_;
    my (%resAtmNames, $i, %atmList, %newRes, $j, $k, $bondIndex);
    my ($resName, @tmp, $node, $type, $atmName, $atmIndex, %COMMON);

    $node = $atomsList->{HEAD};
    $atmName = $node->{data}{ATMNAME};
    $atmIndex = $node->{data}{INDEX};
    if (exists($resAtomList->{$atmIndex})) {
	$resAtmNames{OLD}{$atmName} = $node;
    }
    $node = $node->{next};
    while (! exists($node->{ISHEAD})) {
	$atmName = $node->{data}{ATMNAME};
	$atmIndex = $node->{data}{INDEX};
	if (exists($resAtomList->{$atmIndex})) {
	    $resAtmNames{OLD}{$atmName} = $node;
	}
	$node = $node->{next};
    }
    
    $node = $resList->{HEAD};
    $atmName = $node->{data}{ATMNAME};
    $resName = $node->{data}{RESNAME};
    $resAtmNames{NEW}{$atmName} = $node;
    $node = $node->{next};
    while (! exists($node->{ISHEAD})) {
	$atmName = $node->{data}{ATMNAME};
	$resAtmNames{NEW}{$atmName} = $node;
	$node = $node->{next};
    }

    for $i (keys %{ $resAtmNames{NEW} }) {
	if (! exists($resAtmNames{OLD}{$i})) {
	    $atmList{$i} = $resAtmNames{NEW}{$i};  
	} else { # these are the atoms that are common
	    $COMMON{$i} = $resAtmNames{NEW}{$i}{data};
	    @tmp = grep !/^BOND|MOL|NODE/, keys %{ $resAtmNames{OLD}{$i}{data} };
	    $atmIndex = $resAtmNames{NEW}{$i}{data}{INDEX};
	    for $j (@tmp) {
		$newRes{$atmIndex}{$j} = $resAtmNames{OLD}{$i}{data}{$j};
	    }
	    $newRes{$atmIndex}{RESNAME} = $resName;
	    $newRes{$atmIndex}{RESNUM} = $resNum;
	    $newRes{$atmIndex}{INDEX} = $atmIndex;
	}
    }

    return () if (! %newRes or scalar(keys %newRes ) < 4); # need to have at least 3 atoms in common to continue
    for $i (keys %atmList) {
	&growAtom(\%COMMON, $atmList{$i}, \%newRes);
	$atmIndex = $atmList{$i}{data}{INDEX};
	$newRes{$atmIndex}{RESNUM} = $resNum;
    }

    # add bonds between new residue atoms only
    for $i (keys %atmList) {
	$atmIndex = $atmList{$i}{data}{INDEX};
	@tmp = grep /^BOND/, keys %{ $atmList{$i} };
	for $j (@tmp) {
	    $bondIndex = $atmList{$i}{$j}{data}{INDEX};
	    $j =~ /^BOND(\d+)/;
	    $k = $1;
	    $k++ while (exists($newRes{$atmIndex}{"BOND" . $k}));
	    $newRes{$atmIndex}{"BOND${k}"}{data} = $newRes{$bondIndex};
	}	
    }

    # add bonds between new residue atoms and same old residue atoms
    for $i (keys %{ $resAtmNames{NEW} }) {
	if (exists($resAtmNames{OLD}{$i})) {
	    $atmIndex = $resAtmNames{NEW}{$i}{data}{INDEX};

	    # add the bonds to the unaffected atoms
	    @tmp = grep /^BOND/, keys %{ $resAtmNames{OLD}{$i} };
	    for $j (@tmp) {
		next if (! keys %{ $resAtmNames{OLD}{$i}{$j}{data} });
		$j =~ /^BOND(\d+)/;
		$k = $1;
		$k++ while (exists($newRes{$atmIndex}{"BOND" . $k}));
		$newRes{$atmIndex}{"BOND${k}"} = $resAtmNames{OLD}{$i}{$j};
		$k = 0;
		$k++ while (exists($resAtmNames{OLD}{$i}{$j}{data}{NODE}{"BOND${k}"}));
		$resAtmNames{OLD}{$i}{$j}{data}{NODE}{"BOND${k}"}{data} = $newRes{$atmIndex};
	    }

	    @tmp = grep /^BOND/, keys %{ $resAtmNames{NEW}{$i} };
	    for $j (@tmp) {
		$bondIndex = $resAtmNames{NEW}{$i}{$j}{data}{INDEX};
		$j =~ /^BOND(\d+)/;
		$k = $1;
		$k++ while (exists($newRes{$atmIndex}{"BOND" . $k}));
		$newRes{$atmIndex}{"BOND${k}"}{data} = $newRes{$bondIndex};
	    }
	}
    }

    return \%newRes;
}

sub growAtom {
    my ($oldAtoms, $atom, $newAtoms) = @_;
    my ($dim, $v1, $v2, $v3, $PMatrix, @atmList, $inverseMatrix, $CATMS); 
    my ($atomIndex, $basisMatrix, $coordVector, @tmp, $i, $j, $C);

    @atmList = grep !/H/, keys %{ $oldAtoms };
    $CATMS = {0 => $oldAtoms->{ $atmList[0] }, 1 => $oldAtoms->{ $atmList[1] }, 2 => $oldAtoms->{ $atmList[2] }};
    $C = CoM($CATMS);
    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$v1->{$dim} = $oldAtoms->{ $atmList[0] }{$dim} - $oldAtoms->{ $atmList[1] }{$dim} - $C->{$dim};
	$v2->{$dim} = $oldAtoms->{ $atmList[1] }{$dim} - $oldAtoms->{ $atmList[2] }{$dim} - $C->{$dim};
    }
    $v3 = CrossProduct($v1, $v2);
    $basisMatrix = Math::MatrixReal->new_from_rows([
						    [$v1->{XCOORD}, $v1->{YCOORD}, $v1->{ZCOORD}],
						    [$v2->{XCOORD}, $v2->{YCOORD}, $v2->{ZCOORD}],
						    [$v3->{XCOORD}, $v3->{YCOORD}, $v3->{ZCOORD}]
						    ]);
    $coordVector = Math::MatrixReal->new_from_rows([
						    [$atom->{data}{XCOORD},
						     $atom->{data}{YCOORD},
						     $atom->{data}{ZCOORD}]
						    ]);
    $inverseMatrix = $basisMatrix->inverse;
    $PMatrix = $coordVector * $inverseMatrix;

    #get new coords
    for $i (0 .. 2) {
	for $j (keys %{ $newAtoms }) {
	    if ($newAtoms->{$j}{ATMNAME} eq $oldAtoms->{ $atmList[$i] }{ATMNAME}) {
		$atmList[$i] = $j;
		last;
	    }
	}
    }

    $CATMS = {0 => $newAtoms->{ $atmList[0] }, 1 => $newAtoms->{ $atmList[1] }, 2 => $newAtoms->{ $atmList[2] }};
    $C = CoM($CATMS);
    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$v1->{$dim} = $newAtoms->{ $atmList[0] }{$dim} - $newAtoms->{ $atmList[1] }{$dim} - $C->{$dim};
	$v2->{$dim} = $newAtoms->{ $atmList[1] }{$dim} - $newAtoms->{ $atmList[2] }{$dim} - $C->{$dim};
    }
    $v3 = CrossProduct($v1, $v2);
    $basisMatrix = Math::MatrixReal->new_from_rows([
						    [$v1->{XCOORD}, $v1->{YCOORD}, $v1->{ZCOORD}],
						    [$v2->{XCOORD}, $v2->{YCOORD}, $v2->{ZCOORD}],
						    [$v3->{XCOORD}, $v3->{YCOORD}, $v3->{ZCOORD}]
						    ]);
    $coordVector = $PMatrix * $basisMatrix;

    $atomIndex = $atom->{data}{INDEX};
    @tmp = grep !/^BOND|MOL|NODE/, keys %{ $atom->{data} };
    for $i (@tmp) {
	$newAtoms->{$atomIndex}{$i} = $atom->{data}{$i};
    }
    $newAtoms->{$atomIndex}{XCOORD} = $coordVector->element(1,1) - $C->{XCOORD};
    $newAtoms->{$atomIndex}{YCOORD} = $coordVector->element(1,2) - $C->{YCOORD};
    $newAtoms->{$atomIndex}{ZCOORD} = $coordVector->element(1,3) - $C->{ZCOORD};
}


sub findNode {
    my ($atomList, $atmName) = @_;
    my ($node);
    
    $node = $atomList->{HEAD};
    return $node if ($node->{data}{ATMNAME} eq $atmName);
    $node = $node->{next};
    while (! exists($node->{ISHEAD})) {
	return $node if ($node->{data}{ATMNAME} eq $atmName);
	$node = $node->{next};
    }
    return ();
}

sub InsertResidue {
    my ($atoms, $bonds, $resAtoms, $resBonds, $offset, $startAtm) = @_;
    my ($i, $atmIndex, $atmOffset, $j, @atmList);

    for $i (keys %{ $atoms }) {
	if ($i > $startAtm) {
	    $atmIndex = $i + $offset;
	    $atoms->{$i}{NEWINDEX} = $i + $offset;
	}
    }
    for $i (keys %{ $bonds }) {
	for $j (0 .. $#{ $bonds->{$i} }) {
	    if ($bonds->{$i}[$j] > $startAtm) {
		$bonds->{$i}[$j] = $atoms->{ $bonds->{$i}[$j] }{NEWINDEX};
	    }
	}
    }
    for $i (keys %{ $atoms }) {
	if ($i > $startAtm) {
	    $atmIndex = $i + $offset;
	    $atoms->{$i}{NEWINDEX} = $atmIndex;
	    $atoms->{$atmIndex} = $atoms->{$i};
	    #delete $atoms->{$i};
	    #delete $bonds->{$i};
	}
    }

    @atmList = sort numerically keys %{ $resAtoms };
    for $i (@atmList) {
	$atmOffset = $i - $atmList[0];
	$atmIndex = $atmOffset + $startAtm;
	$atoms->{$atmIndex} = $resAtoms->{$i};
	$bonds->{$atmIndex} = ();
	for $j (@{ $bonds->{$i} }) {
	    push @{ $bonds->{$atmIndex} }, ($j - $atmList[0] + $startAtm);
	}
    }
}

sub init {
    my (%OPTS, $mutateInfo);
    getopt('bsrm',\%OPTS);
    
    for ("b", "r") {
	die "usage: $0 -b bgf file -r residue bgf file -m (mutation residue(s) default all) -s (savename)\n"
	    if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    for ("b", "r") {
	FileTester($OPTS{$_});
    }
    ($bgfFile,$saveFile,$resFile,$mutateInfo) = ($OPTS{b},$OPTS{s},$OPTS{r},$OPTS{m});
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_" . basename($resFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_mutate.bgf";
    }
    $mutateInfo = "*" if (! defined($mutateInfo));
    if ($mutateInfo =~ /\s+/) {
        @{ $selection } = split /\s+/, $mutateInfo;
    } else {
        $selection->[0] = $mutateInfo;
    }
    print "Done\n";
}
