#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/home/yjn1818/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::FileFormats qw(GetPDBFileInfo);
use Packages::General qw(FileTester);

sub createPQR;
sub init;
sub ReadJagOut;
sub addRadii;
sub numerically { ($a<=>$b); }

my ($pdbFile, $jagOutFile, $saveFile);
my ($ATOMS, $BONDS, $JAGDATA);

$|++;
&init;
print "Parsing PDBFile $pdbFile...";
($ATOMS, $BONDS) = GetPDBFileInfo($pdbFile);
print "Done\nParsing Jaguar output file $jagOutFile...";
$JAGDATA = ReadJagOut($jagOutFile);
print "Done\nAdding Radii Information to PDB file...";
&addRadii($ATOMS, $JAGDATA);
print "Done\nCreating PQR file $saveFile...";
createPQR($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub createPQR {
    my ($atoms, $bonds, $saveName) = @_;
    my ($atmName, $resName, $i, $bondList, $j);
    my ($fmt) = '%-6s%5d %-5s%3s %5d %11.3f%8.3f%8.3f%8.3f%8.3f';
    open PDBFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    print PDBFILE "HEADER    PROTEIN\nCOMPND    $saveName\nAUTHOR    GENERATED BY $0\n";
    for $i (sort numerically keys %{ $atoms }) {
        $atmName = $atoms->{$i}{"ATMNAME"};
        $resName = $atoms->{$i}{"RESNAME"};
        $atmName = substr($atmName, 0,3) . substr($atmName,-1,1) if (length($atmName) > 4);
        $resName = substr($resName, 0,2) . substr($resName,-1,1) if (length($resName) > 3);
        printf PDBFILE "$fmt\n", "ATOM", $i, $atmName, $resName, $atoms->{$i}{RESNUM},
        $atoms->{$i}{XCOORD}, $atoms->{$i}{YCOORD}, $atoms->{$i}{ZCOORD}, $atoms->{$i}{CHARGE}, 
	$atoms->{$i}{RADII};
        $bondList .= sprintf("%-6s%5d", "CONECT", $i);
        for $j (@{ $bonds->{$i} }) {
            $bondList .= sprintf("%5d",$j);
        }
        $bondList .= "\n";
    }

    print PDBFILE "CONECT\n$bondList" if (keys %{ $bonds });
    close PDBFILE;
}

sub addRadii {
    my ($atoms, $jagData) = @_;
    my ($i, $resNum, $atmName);
    
    for $i (keys %{ $atoms }) {
	$resNum = $atoms->{$i}{RESNUM};
	$atmName = $atoms->{$i}{ATMNAME};
	if (exists($jagData->{$i}{RADII})) {
	    $atoms->{$i}{RADII} = $jagData->{$i}{RADII};
	    $atoms->{$i}{CHARGE} = $jagData->{$i}{CHARGE};
	}
    }
}

sub ReadJagOut(@) {
    my ($jag_out) = $_[0];
    my ($inline, $is_radii, $hash_key, $is_geometry, $is_charge);
    my (@atom_id, @atom_charge, $counter, $is_valid, $curr_id);
    my (%Jaguar_Data, $resCounter);

    $is_radii = $is_geometry = $is_charge =  $is_valid =  0;
    $resCounter = 1;
    open JAGOUT, $jag_out or die "Cannot read from $jag_out: $!\n";
    while (<JAGOUT>) {
        chomp;
        $inline = $_;
        if ($inline =~ /final geometry/) {
            $is_geometry = 1;
        } elsif ($inline =~ /principal moments of inertia/) {
            $is_geometry = 0;
        }elsif ($is_geometry) {
            if ($inline =~ /([A-Z]{1}[a-z]?)(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
                $hash_key = $2;

                $Jaguar_Data{$hash_key}{"XCOORD"} = $3;
                $Jaguar_Data{$hash_key}{"YCOORD"} = $4;
                $Jaguar_Data{$hash_key}{"ZCOORD"} = $5;
		$Jaguar_Data{$hash_key}{"ATMNAME"} = "${1}${2}";
                $is_valid = 1;
            }
        } elsif ($inline =~ /Atomic charges from/) {
            $is_charge = 1;
        } elsif ($inline =~ /sum of atomic charges/) {
            $is_charge = 0;
        }elsif ($is_charge && $is_valid) {
            if ($inline =~ /Atom\s+(.+)$/) {
                @atom_id = split /\s+/, $1;
            }elsif ($inline =~ /Charge\s+(.+)$/) {
                @atom_charge = split /\s+/, $1;
                for $counter (0 .. $#atom_charge) {
                    $atom_id[$counter] =~ /(\d+)/;
                    $curr_id = $1;

                    if ($Jaguar_Data{$curr_id}) {
#                       print "Wrote Charge\n";
                        $Jaguar_Data{$curr_id}{"CHARGE"} = $atom_charge[$counter];
                        $is_valid = 1;
                    } else {
                        print "Charge found for Non Existant Atom: $atom_id[$counter]\n";
                    }
                }
            }
        } elsif ($inline =~ /van der Waals radii/) {
            $is_radii = 1;
            $counter = 0;
        } elsif ($inline =~ /Number of Lewis structures found/i) {
            $is_radii = 0;
        } elsif ($is_radii and $inline =~ /^\s*[A-Z]{1}[a-z]?(\d+)\s+(\d+\.\d+)/) {
            $Jaguar_Data{$1}{"RADII"} = $2;
        }
    }
    close JAGOUT;

    die "Error reading Jaguar Output $jag_out\n" if (! $is_valid);

    return (\%Jaguar_Data);
}

sub init {
    my (%OPTS);
    getopt('pjs',\%OPTS);
    
    for ("p", "j") {
	die "usage: $0 -p pdbfile -j jaguar ouput file -s (pqr file)\n" if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($pdbFile, $jagOutFile, $saveFile) = ($OPTS{p}, $OPTS{j}, $OPTS{s});
    FileTester($pdbFile);
    FileTester($jagOutFile);
    
    if (! defined($saveFile)) {
	$saveFile = basename($pdbFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= ".pqr";
    }
    print "Done\n";
}