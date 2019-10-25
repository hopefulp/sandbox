#!/usr/bin/perl
# written by Joonho Park
# read .msi and write .lxyz .pos(vasp) as cartesian coordinate


if($#ARGV < 0){
    print "Error: input msi input file";
    exit(0);
}

$fmsi=$ARGV[0];
$nspecies=$#ARGV;
@atoms=();
for($i=0;$i<$nspecies;$i++){
   $atoms[$i]=$ARGV[$i+1];
} 

@fmsi=split(/\./,$fmsi);
$fmsi=$fmsi[0].".msi";

print @atoms,"\n";
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

for($i=0;$i<$nspecies;$i++){
    $natom_species=0;
}

@alattice=();
@blattice=();
@clattice=();


$find="NO";
$flag="OFF";
@bas_line=();
### read three atoms's coordinates
#$i=0;
while($line=<IN>){
#    $i++;
    #print $line;
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
	for($i=0;$i<=$#atoms;$i++){
	    $check_species="NO";
	    if($atoms[$i] eq $coord[$natom][3]){
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
	}
#	print $coord[$natom][0],$coord[$natom][1],$coord[$natom][2],$coord[$natom][3],"\n";
	$natom++; 
    }	
}

print "natom = ",$natom,"\n";
print OUTvasp "    1.0000\n";
for($i=0;$i<3;$i++){
    print "\t$alattice[$i]";
    print OUTvasp "\t$alattice[$i]";
}   print OUTvasp "\n"; print "\n";
for($i=0;$i<3;$i++){
    print "\t$blattice[$i]";
    print OUTvasp "\t$blattice[$i]";
}   print OUTvasp "\n"; print "\n";
for($i=0;$i<3;$i++){
    print "\t$clattice[$i]";
    print OUTvasp "\t$clattice[$i]";
}   print OUTvasp "\n"; print "\n";

for($i=0;$i<=$#atoms;$i++){
    print "  $atoms[$i] ";
    print OUTvasp "  $atoms[$i] ";
}   print OUTvasp "\n";  print "\n";

for($i=0;$i<$nspecies;$i++){
    print "  $natom_species[$i]";
    print OUTvasp  "  $natom_species[$i]";
}   print OUTvasp   "\n"; print "\n";

print OUTvasp "Cartesian\n";
for($j=0;$j<=$#atoms;$j++){
    for($i=0;$i<$natom;$i++){
	if($coord[$i][3] eq $atoms[$j]){
	    for($k=0;$k<3;$k++){
		print "\t$coord[$i][$k]";
		print OUTvasp "\t$coord[$i][$k]";
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

