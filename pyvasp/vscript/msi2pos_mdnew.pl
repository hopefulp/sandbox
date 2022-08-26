#!/usr/bin/perl
# written by Joonho Park
# read .msi and write .lxyz .pos(vasp) as cartesian coordinate
# usage: ./msi2pos_mdnew.pl cpo27zn.msi Zn O C H	# atoms are in the order of POTCAR
# usage: ./msi2pos_md.pl  cpo27zn.msi 0 12 6 0 Zn O C H # MD_false is for number of atoms which will be fixed in MD
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

$nspecies=$#ARGV;
@atom_kinds=();
%atom_kinds=();
for($i=0;$i<=$nspecies;$i++){
    if($ARGV[$i] =~ /=/){
    	@list=split(/=/,$ARGV[$i]);
        push(@atom_kinds,$list[0]);
        $atom_kinds{$list[0]}=$list[1];
    }else{
        push(@atom_kinds,$ARGV[$i]);
        $atom_kinds{$ARGV[$i]}=0;
    }
}

### T or False
$dyn = 'T';

$atom_kinds = join("   ",@atom_kinds);
print "   ",$atom_kinds,"\n";

for($i=0;$i<=$#atom_kinds;$i++){
    print "MD atom_kinds are $atom_kinds[$i] $atom_kinds{$atom_kinds[$i]} molecules\n";
}
#exit;

$flxyz=$fmsi[0].".lxyz";
$fvasp="POSCAR.".$fmsi[0];
open(IN,"<$fmsi");
open(OUTlxyz,">$lxyz");
open(OUTvasp,">$fvasp");

print OUTvasp @atom_kinds,"\n";
print "Output to ", $fvasp, " ", $flxyz, "\n";

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
@max_kind=qw( 0 0 0 0);
@i_kind=qw( 0 0 0 0 );
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
#	print $natom_kind[$natom];
	$check_species="NO";
	for($i=0;$i<=$#atom_kinds;$i++){
	    if($atom_kinds[$i] eq $natom_kind[$natom]){
		$max_kind[$i]++;
	#	$natom_species[$i]++; 
		$check_species="YES";
		last;
	    }
        }
	if($check_species eq "NO"){
	   print "ERROR: $natom_kind[$natom]  not designated is detected\n";
	   exit;
	}
#	print $field[5];
    }
    if($field[$xyz_label] eq $key_coord){
	for($i=0;$i<3;$i++){
	    $coord[$natom][$i]=$field[$xyz_label+1+$i];
#	    if($natom<$nmd_false){
#		$coord[$natom][$i+3]="F";
#	    }else{
#		$coord[$natom][$i+3]="T";
#	    }
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

#$natom_species = join("   ",@natom_species);
print "   ",$atom_kinds,"\n";
print OUTvasp "   ",$atom_kinds,"\n";

for($i=0;$i<=$#atom_kinds;$i++){
    print "  $max_kind[$i]";
    print OUTvasp "  $max_kind[$i]";
}   print "\n"; print OUTvasp "\n";

print OUTvasp "Selective Dynamics\n" ;
print OUTvasp "Cartesian\n";

if($dyn eq 'T'){
    $str = "  T T T";
}else{
    $str = "  F F F";
}

for($j=0;$j<=$#atom_kinds;$j++){
    for($i=0;$i<$natom;$i++){
	    if($natom_kind[$i] eq $atom_kinds[$j]){
            $i_kind[$j]++;
#           print $max_kind[$j],"\t",$atom_kinds{$atom_kinds[$j]}."\n";
            for($k=0;$k<3;$k++){
                printf "  %10.5f", $coord[$i][$k];
                printf OUTvasp "  %10.5f", $coord[$i][$k];
            }
            if($i_kind[$j] <= ($max_kind[$j]-$atom_kinds{$atom_kinds[$j]})){
                print           "$str\n";
                print OUTvasp   "$str\n";
            }else{
                print OUTvasp "$str\n";
                print OUTvasp "$str\n";
            }
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

