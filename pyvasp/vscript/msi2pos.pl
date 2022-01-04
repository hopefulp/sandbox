#!/usr/bin/perl
# written by Joonho Park
# read .msi and write .lxyz POSCAR_name(vasp) as cartesian coordinate
# usage: msi2pos.pl a.msi Mn O C H
# poscar will be written in the order of argument
# if you miss an atom in the a.msi, it will report a certain atom is missing.

if($#ARGV < 0){
    print "Error: input msi input file";
    exit(0);
}

$fmsi=shift(@ARGV);
$nspecies=$#ARGV+1;
@atoms=();
for($i=0;$i<$nspecies;$i++){
   $atoms[$i]=$ARGV[$i];
} 

@fmsi=split(/\./,$fmsi);
if($fmsi[$#fmsi] eq "msi"){
    if($#fmsi == 1){
    	$fname=$fmsi[0];
    }else{
	$fname="";
	for($i=0;$i<$#fmsi;$i++){
	    if($i==0){  $fname.="$fmsi[$i]";
	    }else{ 	$fname.=".$fmsi[$i]";}
	}
    }
}else{
    print "It is not msi file\n";
    exit(1);
}

$atoms = join("   ",@atoms);
print "   ",$atoms,"\n";
#exit;

$flxyz=$fname.".lxyz";
$fvasp="POSCAR.".$fname;
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


$find="NO";
$flag="OFF";
@bas_line=();
### read three atoms's coordinates
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
	$coord[$natom][3]=$field[5];
	$check_species="NO";
	$detected_atom=$coord[$natom][3];
	for($i=0;$i<=$#atoms;$i++){
	    if($atoms[$i] eq $detected_atom){
		$natom_species[$i]++; 
		$check_species="YES";
		last;
	    }
        }
	if($check_species eq "NO"){
	   print "ERROR: atom not designated is detected \"$detected_atom\"\n";
	   exit;
	}
#	print $field[5];
    }
    if($field[$xyz_label] eq $key_coord){
	for($i=0;$i<3;$i++){
	    $coord[$natom][$i]=$field[$xyz_label+1+$i];
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

print OUTvasp "Cartesian\n";
for($j=0;$j<=$#atoms;$j++){
    for($i=0;$i<$natom;$i++){
	if($coord[$i][3] eq $atoms[$j]){
	    for($k=0;$k<3;$k++){
		printf "\t%10.5f", $coord[$i][$k];
		printf OUTvasp "\t%10.5f", $coord[$i][$k];
	    }   print "  $coord[$i][3]\n"; print OUTvasp "\n";
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

