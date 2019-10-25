#!/usr/bin/perl5.8.8 -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use Packages::MolData;
use Getopt::Std qw(getopt);

sub init;
sub showUsage;
sub getSolvent;
sub getCellOpts;
sub getSolvOpts;
sub createSolventBox;
sub embedMols;
sub removeSolv;
sub getSoluData;

my ($sfile, $ffType);
my ($solute, $solvent, $readFunc, $ffields, $cellOpts, $system, $tot, $totres);

$|++;
&init;
print "Gettting data from solute $sfile->{type} file $sfile->{name}...";
$solute->read($sfile->{name}, $sfile->{type});
print "Done\nCreating solvent box...";
&createSolventBox($solvent, $solute, $cellOpts);
$solute->write($sfile->{save}, "bgf");
$tot = getSoluData($solute);
undef($solute);
print "Done\nEmbedding solute and removing bad contacts...";
&embedMols($sfile, "_solv_trim.bgf", $ffields);
print "Done\n";#Loading embeded system...";
#$system =  Packages::MolData->new();
#$system->read($sfile->{save}, "bgf");
#print "Done\n";
#$tot = &removeSolv($system, $cellOpts, $tot);
$tot = removeSolv($sfile->{save}, $cellOpts, $tot);
print "Added $tot->{atoms} solvent atoms ($tot->{molecule} molecules)...Done\n"; 
#print "Creating $sfile->{type} file $sfile->{save}...";
#$system->write($sfile->{save}, $sfile->{type});
#print "Done\n";

sub getSoluData {
    my ($sol) = $_[0];
    my ($tot, @tmp);

    %{ $tot } = %{ $sol->count };
    @tmp = sort {$a<=>$b} keys %{ $sol->shash->{resid} };
    $tot->{resid} = pop @tmp;
    return $tot;
}

sub removeSolv {
    my ($bgffile, $cell, $soluteCount) = @_;
    my ($tot, $count, $solvCount, $countStr, $getAtmStr, $start, $i);

    $start = $soluteCount->{resid};
    print "Computing stats...";
    $countStr = "/home/yjn1818/scripts/countAtoms.pl -f $bgffile -t bgf -s \"resnum>${start}\" -m 1";
    open COUNT, "$countStr |" or die "ERROR while executing $countStr\n";
    while (<COUNT>) {
	if ($_ =~ /^Found (\d+) (atoms|molecule)/) {
	    $count->{$2} = $1;
	}
    }
    close COUNT;
    $tot = $count->{molecule}; # total number of solvent molecules
    if (! exists($cell->{total}) or $tot < $cell->{total}) {
        for $i ("atoms", "molecule") {
            $solvCount->{$i} = $count->{$i} - $soluteCount->{$i};
        }
        return $solvCount;
    }

    $tot -= $cell->{total};
    print "Done\nRemoving $tot solvent molecules to achieve $cell->{total}...";
    $getAtmStr = "$Bin/removeMols.pl -b $bgffile -s $bgffile -a \"resnum>${start}\" -m $tot > _out.dat";
    die "\n" if (system($getAtmStr));
    open COUNT, "$countStr |" or die "ERROR while executing $countStr\n";
    while (<COUNT>) {
        if ($_ =~ /^Found (\d+) (atoms|molecule)/) {
            $count->{$2} = $1;
        }
    }
    close COUNT;

    print "Done\n";
    return $count;
}
sub removeSolv_old {
    my ($structure, $cell, $count) = @_;
    my ($solvCount, $i, $tot, $mol, $molId, $start, @molIDs);
    
    $tot = $structure->count->{molecule} - $count->{molecule}; # total number of solvent molecules
    if (! exists($cell->{total}) or $tot < $cell->{total}) {
	for $i ("atoms", "molecule") {
	    $solvCount->{$i} = $structure->count->{$i} - $count->{$i};
	}
        return $solvCount;
    }

    $start = $count->{molecule};
    $tot -= $cell->{total};
    for $i ($start .. $structure->count->{molecule}) {
	push @molIDs, $i;
    }
    $i = 1;
    while ($i <= $tot) {
	$molId = int(rand($#molIDs));
	print "Deleting molecule $molIDs[$molId] ($i of $tot)...\r";
	$mol = $structure->molecule($molIDs[$molId]);
	splice @molIDs, $molId, 1;
	$structure->deleteMol($mol);
	$i++;
    }
    print "Deleted $tot random solvent molecules to get to $cell->{total}..\n";
    for $i ("atoms", "molecule") {
        $solvCount->{$i} = $structure->count->{$i} - $count->{$i};
    }
    return $solvCount;
}

sub embedMols {
    my ($soluFile, $solvFile, $ffields) = @_;
    my ($embed);

    $embed = "/home/yjn1818/scripts/embedMolecule.pl -m $solvFile -s $soluFile->{name} -f \"$ffields\" -c none -w $soluFile->{save}";
    die "\n" if (system("$embed > _out"));
    system("rm -fr _out _solv_1cell.bgf _solv_replicate.bgf _solv_trim.bgf");
}

sub createSolventBox {
    my ($solv, $solu, $cell) = @_;
    my ($i, $replicate, $blen, $trim, $box, $bmin, $smin, $offset, $map, $solvBox, $cellScale);
    my ($tmp, $replicateDim);

    $map = ({ "a" => "x", "b" => "y", "c" => "z"});
    if (defined($cell->{density})) { #compress/expand the solvent cell to the new density. assume 1 g/cm3
	$cellScale = 1/($cell->{density}**(1/3)); #
	for $i ("x", "y", "z") {
	    $solv->stressCell($i, $cellScale);
	}
    }
    $smin = $solv->getExtrema("min");
    %{ $solvBox } = %{ $solv->vdwbox };
    %{ $solvBox } = %{ $solv->cell } if ($solv->cell->{valid});
    if($solu->cell->{valid}) {
    %{ $box }  = %{ $solu->cell };
    for ("a", "b", "c") {
        $tmp = ();
        $tmp = $box->{$_};
        delete $box->{$_};
        $box->{$_}{max} = $tmp;
        $box->{$_}{min} = 0;
    }
    } else {
    %{ $box } = %{ $solu->vdwbox };
	#this is for the case where the solute is small
	#inflate by 1 water shell to ensure that at least
	#some waters remain
	for $i ("a", "b", "c") {
	    $box->{$i}{max} += 7;
	    $box->{$i}{min} -= 7;
	}
    }
    $replicate = "/home/yjn1818/scripts/replicate.pl -b ./_solv_1cell.bgf -d \"";
    $replicateDim =""; 
    for $i ("a", "b", "c") {
	$box->{$i}{max} += $cell->{cell}{$i}{max} 
		if (exists($cell->{cell}) and exists($cell->{cell}{$i}) and exists($cell->{cell}{$i}{max}));
        $box->{$i}{min} -= $cell->{cell}{$i}{min} 
		if (exists($cell->{cell}) and exists($cell->{cell}{$i}) and exists($cell->{cell}{$i}{min}));
	$box->{$i}{len} = $box->{$i}{max} - $box->{$i}{min};
	$bmin->{$i} = $box->{$i}{min};
	$blen .= "$box->{$i}{len} ";
	$replicateDim .= sprintf("%.0f ", (($box->{$i}{len}/$solvBox->{$i}))+1);
    }
    $replicate .= $replicateDim;
	$replicate .= "\" -s _solv_replicate.bgf";
    #move the solvent box to the solute minima
    for $i ("a", "b", "c") {
	$offset->{ $map->{$i }} = $bmin->{$i} - $smin->{$i};
    }
    $solv->moveMol("all", $offset);
    $solv->write("_solv_1cell.bgf", "bgf");
    #replicate the cell by the replication vector calculated above
    if($replicateDim !~ /1 1 1/) {
	die "ERROR while executing \"$replicate\"\n" if (system("${replicate} >& _out.dat"));
	} else {
	system("cp ./_solv_1cell.bgf _solv_replicate.bgf");
	}
    # remove all molecules outside the solute (inflated) cell
	$trim = "/home/yjn1818/scripts/trimCell.pl -b _solv_replicate.bgf -c \"$blen\" -s _solv_trim.bgf -m 1 -o 2";
    die "ERROR while executing \"$trim\"\n" if (system("${trim} >& _out.dat"));
    undef($solv);
}

sub init {
    my (%OPTS, $ffStr, $solvOptsStr, $solvTypeStr, $periodStr, $solvOpts);

    getopt('ifntws', \%OPTS);
    for ("f", "i","n") {
	die &showUsage . "\n" if (! exists($OPTS{$_}));
    }

    print "Initialzing...";
    ($sfile->{name}, $sfile->{type}, $ffields, $solvOptsStr, $solvTypeStr, $sfile->{save}) = 
	($OPTS{i}, $OPTS{t}, $OPTS{f}, $OPTS{n}, $OPTS{w}, $OPTS{s});
    $solute =  Packages::MolData->new();
    $solute->testFile($sfile->{name});
    if (! defined($sfile->{type})) {
	$sfile->{type} = "bgf";
	if ($sfile->{type} =~ /\.(\w+)$/) {
	    $sfile->{type} = lc $1;
	}
    }
    $cellOpts = getCellOpts($solvOptsStr);
    $solvOpts = getSolvOpts($solvTypeStr);
    $sfile->{save} = $solute->getFileName($sfile->{name}) if (! defined($sfile->{save}));
    print "Done\n";
    $solvent = Packages::MolData->new();
    print "Gettting data from solvent $solvOpts->{type} file $solvOpts->{file}...";
    $solvent->read($solvOpts->{file}, $solvOpts->{type});
    print "Done\n";
}

sub getCellOpts {
    my ($solventStr) = $_[0];
    my ($SOLVENT, $dim, $i, $operator, $val, $map);

    $map = ({ "x" => "a", "y" => "b", "z" => "c" });
    if ($solventStr =~ /total:\s*(\d+)/i) {
	$SOLVENT->{total} = $1;
    } elsif ($solventStr =~ /density:\s*(\d+\.?\d*)/i) {
	$SOLVENT->{density} = $1;
    }

    while ($solventStr =~ /(x|y|z|x.|y.|z.|x..|y..|z..):\s*(\+|\-|\+\/?\-)\s*(\d+\.?\d*)/gi) {
	($dim, $operator, $val) = (lc $1, "+", $2);
	($operator, $val) = ($2, $3) if ($3);
	while ($operator =~ /(\+|\-)/g) {
	    $i = "max" if ($1 eq "+");
	    $i = "min" if ($1 eq "-");
	    while ($dim =~ /(x|y|z)/g) {
		$SOLVENT->{cell}{ $map->{$1} }{$i} = $val;
	    }
	}
    }
    die "ERROR: Unindentified options obtained for solvent options: \"$solventStr\"\n"
	if (! $SOLVENT);
    return $SOLVENT;
}

sub getSolvOpts {
    my ($solventStr) = lc $_[0];
    my ($SOLVENT);

    if (! defined($solventStr)) {
        $SOLVENT->{file} = "/home/yjn1818/scripts/dat/WAT/f3c_box.bgf";
        $SOLVENT->{type} = "bgf";
    }elsif (-e $solventStr and -r $solventStr and -T $solventStr) {
	$SOLVENT->{file} = $solventStr;
	$SOLVENT->{type} = "bgf";
	if ($solventStr =~ /\.(\w+)$/) {
	    $SOLVENT->{type} = lc $1;
	}
    } elsif ($solventStr =~ /(tip3|tip3_charmm|f3c|meso|spc|tip4)/) {
	$SOLVENT->{file} = "/home/yjn1818/scripts/dat/WAT/${1}_box.bgf";
	$SOLVENT->{type} = "bgf";
    }

    if (! $SOLVENT) {
        $SOLVENT->{file} = "/home/yjn1818/scripts/dat/WAT/f3c_box.bgf";
        $SOLVENT->{type} = "bgf";
    }

    return $SOLVENT;
}

sub showUsage {
    my $usage = <<DATA;
usage: $0 -i input structure -f forcefield(s) -n solvent options -t (file type) -w (solvent type) -s (savename)
OPTIONS:
	-i input structure: REQUIRED. The following file formats are supported: BGF, PDB, MSI, MOL2
	-f forcefield(s):   REQUIRED. The solute and corresponding solvent atoms must all have vdw 
			    terms defined in this file(s). Can be either a CERIUS2 or MPSIM formatted forcefield.
	-n solvent options: REQUIRED. Specifies the number of solvent molecules to add. Can be either:
			    "density: x.x" - will add solvent atoms to the cell until density x.x is achieved. Will
					     assume that input solvent is equilibrated at denisty of 1 gm/cm3
			    "total: x"     - will add a total of x solvent molecules
			    "[x|y|z]:[+|-] x.x" - will inflate the unit cell by x.x angstoms in specified direction 
					   (+, -, +/- directions) for the specified dimension (any combo of x,y,z). 
					   This option can be used with the previous 2. If none of the other 2 is 
					   specified, will assume density 1 gm/cm3. Multiple entries should be
					   enclosed in quotes. e.g. "x: +/- 10 y: +10 z: -12"
	-t file type:       OPTIONAL. Specifies the formatting of the input file. If not supplied, with default 
			    to either the file suffix or to a BGF file if the file has no suffix
	-w solvent type:    OPTIONAL. The following predetermined (equilibrated) solvents are available. Default F3C.
			    Alternatively, you can provide the location of your solvent file.
			    TIP3: the original Jorgenson TIP3 water model (rigid hydrogens)
			    TIP4: TIP4P with massless pseudo-atom
			    TIP3_CHARMM: TIP3 water model as implemented in CHARMM
			    F3C: F3C water model (no rigid hydrogens)
			    SPC: SPC water model
			    MESO: Water model for Meso-scale simulations (Valeria Molinero)
	-s savename:        OPTIONAL. Will assume \$prefix_mod.\$suffix if not specified.
DATA

    return $usage;
}
