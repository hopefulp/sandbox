package Packages::ManipAtoms;
require Exporter;
use strict;
use Cwd;

our (@ISA, @EXPORT, @EXPORT_OK, $VERSION);

@ISA = qw(Exporter);
$VERSION = "1.00";
@EXPORT = ();
@EXPORT_OK = qw(ImageAtoms ScaleAtoms UnwrapAtoms FindElement CenterSystem GetAtmList 
	     GetSolvent SplitAtomsByMol GetAtmData GetMols MapOrigin);

sub numerically { ($a<=>$b); }

sub MapOrigin {
    my ($box, $com) = @_;
    my ($boxCenter, $offset);

    for (keys %{ $box }) {
        $boxCenter->{$_} = ($box->{$_}{hi} - $box->{$_}{lo})/2;
        $offset->{$_} = $boxCenter->{$_} - $com->{$_ . "COORD"};
        $box->{$_}{hi} -= $offset->{$_};
        $box->{$_}{lo} -= $offset->{$_};
    }
}

sub CenterSystem {
    my ($atoms, $offset, $atomList) = @_;
    my ($i, $dim);

    for $i (keys %{ $atoms }) {
	next if (keys %{ $atomList } && ! exists($atomList->{$i}));
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    $offset->{$dim} = $offset->{substr($dim,0,1)}{lo} if (! exists($offset->{$dim}));
	    next if (! $offset->{$dim});
	    $atoms->{$i}{$dim} -= $offset->{$dim};
	}
    }
}

sub FindElement {
    my ($fftype, $parms, $elements) = @_;
    my ($elementName, $elementNum, $i);

    die "ERROR: $fftype is not a valid force field type!\n" if (! exists($parms->{$fftype}));
    $elementName = $parms->{$fftype}{ATOM};
    $elementNum = 0;

    for $i (keys %{ $elements }) {
        if (uc($elements->{$i}{SYMBOL}) eq uc($elementName)) {
            $elementNum = $i;
            last;
        }
    }

    #die "ERROR: Element $elementName is not a valid element!\n" if (! $elementNum);

    return ($elementName, $elementNum);
}

sub ImageAtoms {
    my ($atoms, $CoM, $box) = @_;
    my ($MOL, $i, $OFFSET, $dim, $index, $CENTER);

    %{ $CENTER } = %{ $CoM };
    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$OFFSET->{$dim} = 0;
	while ($CENTER->{$dim} > $box->{$dim}{hi}) {
	    $OFFSET->{$dim}++;
	    $CENTER->{$dim} -= $box->{$dim}{len};
        }
	while ($CENTER->{$dim} < $box->{$dim}{lo}) {
	    $OFFSET->{$dim}--;
	    $CENTER->{$dim} += $box->{$dim}{len};
	}
    }
    
    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$index = $dim;
	$index =~ s/COORD/INDEX/;
	for $i (keys %{ $atoms }) {
	    $atoms->{$i}{$dim} -= ($OFFSET->{$dim} * $box->{$dim}{len});
	    $atoms->{$i}{$index} = $OFFSET->{$dim};
	    if ($atoms->{$i}{$dim} < $box->{$dim}{lo}) {
		$atoms->{$i}{$index}--;
	    } elsif ($atoms->{$i}{$dim} > $box->{$dim}{hi}) {
		$atoms->{$i}{$index}++;
	    }
	    $atoms->{$i}{IMAGED} = 1;
	}
    }

    return $CENTER;
}

sub ScaleAtoms {
    my ($atoms, $box) = @_;
    my ($i, $coord, $dim, $index, $pos);

    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$index = $dim;
	$index =~ s/COORD/INDEX/;
	for $i (keys %{ $atoms }) {
	    $pos = $atoms->{$i}{$dim};
	    $pos /= $box->{$dim}{len};
	    if ($pos =~ /(\d+)(\.\d+)/) {
		if ($1 > 0) {
		    $atoms->{$i}{$index} = $1 + 1;
		} else {
		    $atoms->{$i}{$index} = $1;
		}
		$atoms->{$i}{$dim} = sprintf("%.6f ", $pos);
	    }
	}
    }
}

sub UnwrapAtoms {
    my ($ATOMS, $BOX, $scaled) = @_;
    my ($atomC, $atom, $coord, $index, $pos);

    for $atomC (keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$atomC} };
	for $coord ("XCOORD", "YCOORD", "ZCOORD") {
	    $index = $coord;
	    $index =~ s/COORD/INDEX/;
	    $index = $atom->{$index};
	    if (! defined($scaled) or $scaled == 1) {
	        $atom->{$coord} *= $BOX->{$coord}{len};
	        $atom->{$coord} += $BOX->{$coord}{lo};
	    }
	    $atom->{$coord} += ($index * $BOX->{$coord}{len});
	}
    }
}

sub GetAtmList {
    my ($select, $ATOMS) = @_;
    my ($i, %LIST, $field, $val, $rec, $operator, $excluded, $tmp);
    my ($count);

    for $field (keys %{ $select }) {
	$tmp = "";
	for $val (keys %{ $select->{$field} }) {
	    $excluded = 0;
	    $operator = "";
	    if ($val eq "*") {
		$operator = ">";
		$val = "0 and \$ATOMS->{\$i}{INDEX} < 999999";
	    } elsif ($val =~ /^\d+/) { #integer
		$operator = "==";
	    } elsif ($val =~ /^\!(\d+)/) {
		$operator = "!=";
		$val = $1;
		$excluded = 1;
	    } elsif ($val =~ /^(>|<|=|\+|\-)(\d+\.*\d*)/) {
		$operator = $1;
		$val = $2;
	    } elsif ($val =~ /^\!(>|<|=)(\d+\.*\d*)/) {
		$excluded = 1;
		if ($1 eq ">") {
		    $operator = "<";
		} elsif ($1 eq "<") {
		    $operator = ">";
		} else {
		    $operator = "!=";
		}
		$val = $2;
	    } elsif ($val =~ /^(\w+)/) {
		$operator = "eq";
		$val = "\"${1}\"";
	    } elsif ($val =~ /^\!(\w+)/) {
		$excluded = 1;
		$operator = "ne";
		$val = "\"${1}\"";
	    }
	    next if (! $operator);
	    if ($tmp) {
		if (! $excluded) {
		    $tmp .= "or \$ATOMS->{\$i}{${field}} $operator $val ";
		} else {
		    $tmp .= "and \$ATOMS->{\$i}{${field}} $operator $val ";
		}
	    } else {
		$tmp = "{${field}} $operator $val ";
	    }
	}
        if ($rec) {
            $rec .= "and (\$ATOMS->{\$i}$tmp) ";
        } else {
            $rec = "$tmp";
        }
    }

    $count = 0;
    for $i (keys %{ $ATOMS }) {
	$excluded = 0;
	if (! eval('$ATOMS->{$i}' . $rec)) {
	    $excluded = 1;
	}
	if (! $excluded) {
	    $count++;
	    $LIST{$i} = $count;
	}
    }
    
    return \%LIST;
}

sub GetSolvent {
    my ($atoms, $solvType) = @_;
    my (%SOLVOPTS, %SOLVATMS, $i);

    %SOLVOPTS = (
			"WATER" => {
					"MOLSIZE" => 3,
					"RESNAME" => "WAT|HOH|RES",
					"FFTYPE"  => "OW|HW|OT|HT|HO|OH|H_|O_3"
				    }
		);

    return () if (! defined($SOLVOPTS{uc($solvType)}));

    for $i (keys %{ $atoms }) {
	if ($atoms->{$i}{MOLSIZE} == $SOLVOPTS{$solvType}{MOLSIZE}) {
	    if ($atoms->{$i}{RESNAME} =~ /$SOLVOPTS{$solvType}{RESNAME}/ or
		$atoms->{$i}{FFTYPE} =~ /$SOLVOPTS{$solvType}{FFTYPE}/) {
		    $atoms->{$i}{IS_SOLVENT} = 1;
		    $SOLVATMS{$i} = 1;
	    }
	}
    }

    return \%SOLVATMS;
}

sub GetMols {
    my ($ATOMS, $BONDS) = @_;
    my ($atomC, @tmp, $rec, $i, $min, $molNum);

    for $atomC (keys %{ $ATOMS }) {
	delete $ATOMS->{$atomC}{MOLECULE};
	delete $ATOMS->{$atomC}{MOLECULEID};
    }
    
    @tmp = sort numerically keys %{ $ATOMS };
    for $atomC (@tmp) {
	$rec = ();
	$rec->{MEMBERS}{$atomC} = 1;
	if (! $BONDS->{$atomC} || $#{ $BONDS->{$atomC} } == -1) { #ions
	    $ATOMS->{$atomC}{MOLECULE} = \%{ $rec };
	    $molNum++;
	} else {
	    $min = $atomC;
	    for $i (@{ $BONDS->{$atomC} })  {
		$rec->{MEMBERS}{$i} = 1;
		$min = $i if ($i < $min);
	    }
	    if ($min < $atomC) { # found a member which has a lower index, so merge this rec with it, and update
		for $i (keys %{ $rec->{MEMBERS} }) {
		    $ATOMS->{$min}{MOLECULE}{MEMBERS}{$i} = 1;
		    next if ($min == $i or $i == $atomC);
		    if (exists($ATOMS->{$i}{MOLECULE}{MEMBERS})) {
			for (keys %{ $ATOMS->{$i}{MOLECULE}{MEMBERS} }) {
			    $ATOMS->{$min}{MOLECULE}{MEMBERS}{$_} = 1;
			}
		    }
		}
		$rec = \%{ $ATOMS->{$min}{MOLECULE} };
	    } else {
		$rec->{INDEX} = $molNum;
		$molNum++;
	    }

	    for $i (keys %{ $rec->{MEMBERS} })  {
		next if ($i == $min);
		$ATOMS->{$i}{MOLECULE} = \%{ $rec };
	    }
	}
    }

    $molNum = 1;
    for $atomC (@tmp) {
	if (! exists($ATOMS->{$atomC}{MOLECULE}{INDEX})) {
	    $ATOMS->{$atomC}{MOLECULE}{INDEX} = $molNum;
	    $molNum++;
	}
	$ATOMS->{$atomC}{MOLECULEID} = $ATOMS->{$atomC}{MOLECULE}{INDEX};
	$ATOMS->{$atomC}{MOLECULE}{SIZE} = scalar(keys %{ $ATOMS->{$atomC}{MOLECULE}{MEMBERS} });
	$ATOMS->{$atomC}{MOLSIZE} = $ATOMS->{$atomC}{MOLECULE}{SIZE};
    }
}

sub SplitAtomsByMol {
    my ($atoms, $selList) = @_;
    my ($i, %molList, $counter, $j, @tmp);

    $selList = $atoms if (! defined($selList));
    @tmp = sort { ($a<=>$b) } keys %{ $selList };

    for $i (@tmp) {
	next if (exists($atoms->{$i}{IS_SPLIT}));
	for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
	    next if (! exists($selList->{$j}));
	    $atoms->{$j}{IS_SPLIT} = 1;
	    $counter = $atoms->{$j}{MOLECULEID};
	    $molList{$counter}{$j} = $atoms->{$i}{MOLECULE}{MEMBERS}{$j};
	}
    }

    for $i (@tmp) {
        delete $atoms->{$i}{IS_SPLIT}
    }

    return \%molList;
}

sub GetAtmData {
    my ($allAtoms, $atomList) = @_;
    my (%ATOMS, $i);
                                                                                                                                      
    for $i (keys %{ $atomList }) {
        $ATOMS{$i} = \%{ $allAtoms->{$i} };
    }
    return \%ATOMS;
}

1;
