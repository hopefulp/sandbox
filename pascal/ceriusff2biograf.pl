#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Packages::CERIUS2 qw(parseCerius2FF);
use Packages::General qw(FileTester LoadElements);
use File::Basename;
use Getopt::Std qw(getopt);

#use Packages::BIOGRAF qw (saveBiografFF);

die "usage: $0 cerius2ff [saveName]\n"
    if (! @ARGV);

sub init;
sub saveBiografFF;
sub matchElement;
sub stringOrder;
sub findAddedH;
sub getType;

my ($FFs, $saveName);
my ($header, $ELEMENTS, $PARMS, $i);
$ELEMENTS = LoadElements;

$|++;
&init;
for $i (@{ $FFs }) {
    print "Parsing Cerius2 FF $i...";
    $PARMS = parseCerius2FF($i, 0, $PARMS);
    print "Done\n";
}
print "Creating Biograf FF $saveName...";
saveBiografFF($PARMS, $saveName, $header);
print "Done\n";

sub init {
    my (%OPTS, $ffStr);

    getopt('fs',\%OPTS);
    die "usage: $0 -f \"cerius2 ff(s)\" -s (savename)\n" if (! exists($OPTS{f}));
    print "Initializing...";
    ($ffStr, $saveName) = ($OPTS{f}, $OPTS{s});
    while ($ffStr =~ /(\S+)/g) {
	push @{ $FFs }, $1 if (-e $1 and -r $1 and -T $1);
    }
    die "ERROR: No valid CERIUS2 forcefields found while searching \"$ffStr\"!\n" if (! defined($FFs));
    if (! defined($saveName)) {
	$saveName = basename($FFs->[0]);
	$saveName =~ s/\.\w+$/_biograf\.par/;
    }

    open HEADER, "/home/yjn1818/scripts/dat/biograf_header.dat" || die "ERROR: Cannot open header file!\n";
    while (<HEADER>) {
	$header .= $_;
    }
    close HEADER;
    print "Done\n";
}

sub saveBiografFF {
    my ($parms, $save, $header) = @_;
    my ($i, $ele, $curr, $addedH, @tmp, $junk); 
    my ($j, $k, $l, $val, $t, $counter);

    @tmp = sort stringOrder keys %{ $parms->{ATOMTYPES} };
    open BIOGRAF, "> $save" || die "ERROR: Cannot create Biograf FF $save: $!\n";
    # FF Options
    print BIOGRAF "$header";
    print BIOGRAF "FFLABEL    ATNO      MASS CHARG HYB BND CPK \#IH \#LP\n";
    for $i (@tmp) {
	$curr = $parms->{ATOMTYPES}{$i};
	$ele = matchElement($curr->{ATOM});
	printf BIOGRAF "%-14s%-4d%-10.4f%-6.1d%-4d%-4d%-4d%-4d%-4d\n", $i, $ele, $curr->{MASS},
	$curr->{CHARGE}, $curr->{NUMBONDS}, 3, 5, $curr->{OTHER}, $curr->{LONEPAIRS};
    }
    print BIOGRAF "*\n";

    #Added H
    print BIOGRAF "ADDED H   HYDROGEN  1IMPLCTH  2IMPLCTH  3IMPLCTH  4IMPLCTH\n";
    $addedH = findAddedH($parms->{BONDS});
    print BIOGRAF "*\nLONEPAIRS\n*\nGASTEIGER          A         B         C        X+\n*\n";
    
    # VDWS
    print BIOGRAF "VDW AT ITY       RNB      DENB     SCALE\n";
    for $i (@tmp) {
	$curr = $parms->{VDW}{$i}{$i}{1};
	printf BIOGRAF "%-6s%4d", $i, getType($curr->{TYPE}, "VDW");
	for $val (@{ $curr->{VALS} }) {
	    printf BIOGRAF "%10.4f", $val;
	}
	print BIOGRAF "\n";
    }
    print BIOGRAF "*\nAUTOTYPE  ELEMENT   HYBRIDIZATION RING_SIZE  REQUIREMENTS  FACTOR\n*\n";
    
    # BONDS
    print BIOGRAF "BONDSTRTCH  TYPE FORC CNST  BND DIST    DE/CUB         E         F  BOND DIP\n";
    for $i (sort stringOrder keys %{ $parms->{BONDS} }) {
	for $j (sort stringOrder keys %{ $parms->{BONDS}{$i} }) {
	    @{ $junk } = keys %{ $parms->{BONDS}{$i}{$j} };
	    $curr = $parms->{BONDS}{$i}{$j}{ shift @{ $junk } };
	    printf BIOGRAF "%-5s-%-5s%5d", $i, $j, getType($curr->{TYPE}, "BONDS");
	    for $val (@{ $curr->{VALS} }) {
		printf BIOGRAF "%10.4f", $val;
	    }
	    print BIOGRAF "\n";
	}
    }
    print BIOGRAF "*\n";

    # ANGLES
    print BIOGRAF "ANGLE             TYPE FORC CNST EQUIL ANG         D         E         F\n";
    for $i (sort stringOrder keys %{ $parms->{ANGLES} }) {
	for $j (sort stringOrder keys %{ $parms->{ANGLES}{$i} }) {
	    for $k (sort stringOrder keys %{ $parms->{ANGLES}{$i}{$j} }) {
		@{ $junk } = keys %{ $parms->{ANGLES}{$i}{$j}{$k} };
		$curr = $parms->{ANGLES}{$i}{$j}{$k}{ shift @{ $junk } };
		printf BIOGRAF "%-5s-%-5s-%-5s%5d", $i, $j, $k, getType($curr->{TYPE}, "ANGLES");
		for $val (@{ $curr->{VALS} }) {
		    printf BIOGRAF "%10.4f", $val;
		}
		print BIOGRAF "\n";
	    }
	}
    }
    print BIOGRAF "*\n";

    # TORSIONS
    $counter = 0;
    print BIOGRAF "TORSION                 CASE   BARRIER    PERIOD CISMIN(1)\n";
    for $i (sort stringOrder keys %{ $parms->{TORSIONS} }) {
	for $j (sort stringOrder keys %{ $parms->{TORSIONS}{$i} }) {
	    for $k (sort stringOrder keys %{ $parms->{TORSIONS}{$i}{$j} }) {
		for $l (sort stringOrder keys %{ $parms->{TORSIONS}{$i}{$j}{$k} }) {
		    $counter++;
		    @{ $junk } = keys %{ $parms->{TORSIONS}{$i}{$j}{$k}{$l} };
		    $curr = $parms->{TORSIONS}{$i}{$j}{$k}{$l}{shift @{ $junk } };
		    $t = 0;
		    while ($t <= ($#{ $curr->{VALS} } - 2)) {
			printf BIOGRAF "%-5s-%-5s-%-5s-%-5s%5d", $i, $j, $k, $l, $counter;
			$curr->{VALS}[$t] *= 4;
			if ($curr->{VALS}[($t + 2)] == 180) {
			    $curr->{VALS}[($t + 2)] = 1;
			} else {
			    $curr->{VALS}[($t + 2)] = -1;
			}
			for $val (0 .. 2) {
			    printf BIOGRAF "%10.4f", $curr->{VALS}[($t + $val)];
			}
			print BIOGRAF "\n";
			$t += 3;
		    }
		}
	    }
	}
    }
    print BIOGRAF "*\n";

    # INVERSIONS
    print BIOGRAF "INVERSION (CENT AT 1ST) TYPE  FRC CNST  EQU ANGL         D         E         F\n";
    for $i (sort stringOrder keys %{ $parms->{INVERSIONS} }) {
	for $j (sort stringOrder keys %{ $parms->{INVERSIONS}{$i} }) {
	    for $k (sort stringOrder keys %{ $parms->{INVERSIONS}{$i}{$j} }) {
		for $l (sort stringOrder keys %{ $parms->{INVERSIONS}{$i}{$j}{$k} }) {
		    @{ $junk } = keys %{ $parms->{INVERSIONS}{$i}{$j}{$k}{$l} };
		    $curr = $parms->{INVERSIONS}{$i}{$j}{$k}{$l}{shift @{ $junk }};
		    printf BIOGRAF "%-5s-%-5s-%-5s-%-5s%5d", $i, $j, $k, $l, getType($curr->{TYPE},"INVERSIONS");
		    for $val (@{ $curr->{VALS} }) {
			printf BIOGRAF "%10.4f", $val;
		    }
		    print BIOGRAF "\n";
		}
	    }
	}
    }
    print BIOGRAF "*\nHBOND       TYPE    -DE HB     RE HB\nEND OF DATA\n";
    close BIOGRAF;
		    
}

sub findAddedH {
    my ($bonds) = $_[0];
    my ($atom1, $atom2, $ele1, $ele2);
    my ($added_h);

    for $atom1 (sort stringOrder keys %{ $bonds }) {
	$ele1 = matchElement($PARMS->{ATOMTYPES}{$atom1}{ATOM});
	for $atom2 (sort stringOrder keys %{ $bonds->{$atom1} }) {
	    $ele2 = matchElement($PARMS->{ATOMTYPES}{$atom2}{ATOM});
	    if ($ele1 == 1) {
		$added_h .= sprintf("%-10s%-10s\n",$atom2, $atom1);
	    } elsif ($ele2 == 1) {
		$added_h .= sprintf("%-10s%-10s\n",$atom1, $atom2);		
	    }
	}
    }

    return $added_h;
}

sub matchElement {
    my ($ele_name) = $_[0];
    my ($ele_num, $i);

    for $i (keys %{ $ELEMENTS }) {
	if ($ELEMENTS->{$i}{SYMBOL} eq $ele_name) {
	    $ele_num = $i;
	    last;
	}
    }

    die "ERROR: No valid element found while search for $ele_name\n"
	if (! defined($ele_num));

    return $ele_num;
}

sub stringOrder {
    ($a cmp $b);
}

sub getType {
    my ($type, $i) = @_;
    my (%TYPES) = (
		   "VDW" => {
		       "LJ_6_12" => 1,
		   },
		   "BONDS" => {
		       "HARMONIC" => 1,
		   },
		   "ANGLES" => {
		       "THETA_HARM" => 21,
		   },
		   "INVERSIONS" => {
		       "IT_JIKL" => 3,
		   },
		   );

    die "ERROR: No type $type in $i found!\n" if (! exists($TYPES{$i}{$type}));
    return $TYPES{$i}{$type};
}
