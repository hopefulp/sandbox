#!/usr/bin/perl
# written by Joonho Park
# read .msi and write .lxyz .pos(vasp) as cartesian coordinate
# usage: ./fmsi2pos_md.pl cpo27zn.msi Zn O C H	# atoms are in the order of POTCAR
# usage: ./fmsi2pos_md.pl  cpo27zn.msi mdf=30 Zn O C H # MD_false is for number of atoms which will be fixed in MD
# note: in the input.msi, md-fixed atoms are to be located in former part and md-movable atoms in later part.

if($#ARGV < 1){
    print "Error: input msi input file\n" if($#ARGV<0);
    print "Error: input atoms symbol array in the order of POTCAR \n";
    exit(0);
}

$fmsi=shift(@ARGV);

@fmsi=split(/\./,$fmsi);
if($fmsi[1] eq ""){
    print "Error: first argument should be filename including .\n";
    exit(0);
}
$fmsi=$fmsi[0].".msi";

print "input file : $fmsi\n";

if($ARGV[0] =~ /=/){
    @nmd=split(/=/,$ARGV[0]);
    $nmd_false=$nmd[1] if ($nmd[0] eq "MDF" || $nmd[0] eq "mdf");
    print "MD false = $nmd_false\n";
    shift(@ARGV);
}

$nspecies=$#ARGV;
@atoms=();
for($i=0;$i<=$nspecies;$i++){
   $atoms[$i]=$ARGV[$i];
} 
$atoms = join("   ",@atoms);
print "   ",$atoms,"\n";
#exit;

$flxyz=$fmsi[0].".lxyz";
$fvasp=$fmsi[0].".pos";
open(IN,"<$fmsi");
open(OUTlxyz,">$lxyz");
open(OUTvasp,">$fvasp");

print OUTvasp @atoms,"\n";


$key_lattice1="A3";
$key_lattice2="B3";
$key_lattice3="C3";
$key_atom="ACL";
$key_coord="XYZ";
$lattice_label=3;
$atom_label=3;
$xyz_label=3;
$natom=0;

@alattice=();
@blattice=();
@clattice=();

@natom_kind=();

$find="NO";
$flag="OFF";
### read file
while($line=<IN>){
#    @field=split(/\s+/,$line);
    @field=split(/[^a-zA-Z0-9\.\-]+/,$line);
#    print @field,"\n";

    if($field[$lattice_label] eq $key_lattice1) {
	for($i=0;$i<3;$i++){
	    push(@alattice,$field[$lattice_label+1+$i]);
	}
    }

    if($field[$lattice_label] eq $key_lattice2) {
	for($i=0;$i<3;$i++){
	    push(@blattice,$field[$lattice_label+$i+1]);
	}
    }
    if($field[$lattice_label] eq $key_lattice3) {
	for($i=0;$i<3;$i++){
	    push(@clattice,$field[$lattice_label+1+$i]);
	}
    }

    if($field[$atom_label] eq $key_atom){
	$natom_kind[$natom]=$field[5];
	$check_species="NO";
	for($i=0;$i<=$#atoms;$i++){
	    if($atoms[$i] eq $natom_kind[$natom]){
		$natom_species[$i]++; 
		$check_species="YES";
		last;
	    }
        }
	if($check_species eq "NO"){
	   print "ERROR: atom not designated is detected\n";
	   exit;
	}
#	print $field[5];
    }
    if($field[$xyz_label] eq $key_coord){
	for($i=0;$i<3;$i++){
	    $coord[$natom][$i]=$field[$xyz_label+1+$i];
	    if($natom<$nmd_false){
		$coord[$natom][$i+3]="F";
	    }else{
		$coord[$natom][$i+3]="T";
	    }
	}
#	print $coord[$natom][0],$coord[$natom][1],$coord[$natom][2],$coord[$natom][3],"\n";
	$natom++; 
    }	
}

print "natom = ",$natom,"\n";
print OUTvasp "    1.0000\n";

$lattice_vector = join("\t",@alattice);
print "\t",$lattice_vector,"\n";
print OUTvasp "\t",$lattice_vector,"\n";
$lattice_vector = join("\t",@blattice);
print "\t",$lattice_vector,"\n";
print OUTvasp "\t",$lattice_vector,"\n";
$lattice_vector = join("\t",@clattice);
print "\t",$lattice_vector,"\n";
print OUTvasp "\t",$lattice_vector,"\n";

print "   ",$atoms,"\n";
print OUTvasp "   ",$atoms,"\n";

$natom_species = join("   ",@natom_species);
print "   ",$natom_species,"\n";
print OUTvasp "   ",$natom_species,"\n";

print OUTvasp "Selective Dynamics\n" if($nmd[0]);
print OUTvasp "Cartesian\n";

for($j=0;$j<=$#atoms;$j++){
    for($i=0;$i<$natom;$i++){
	if($natom_kind[$i] eq $atoms[$j]){
	    for($k=0;$k<3;$k++){
		print "\t$coord[$i][$k]";
		print OUTvasp "\t$coord[$i][$k]";
	    }
	    print OUTvasp "\t$coord[$i][3]  $coord[$i][4]  $coord[$i][5]\n" if($nmd[1]);
	    print "  $coord[$i][3]\n";
 	}
    }
}

#print @bas_line;
#print OUT @bas_line;

#if($find eq "NO"){
#    print "Error: There is not atom $atom in $basfile\n";
#}
close(IN); 
close(OUTvasp);
close(OUTlxyz);

