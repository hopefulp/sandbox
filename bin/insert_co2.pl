#!/usr/bin/perl
# written by Joonho Park
# read  bare mof POSCAR	 and add CO2 extracted from mof-co2 POSCAR and output file

use List::Util qw(sum);

if($#ARGV < 1){
    print "Error: input TWO POSCAR's ";
    exit(0);
}

$fmof=$ARGV[0];
$fco2=$ARGV[1];

if($#ARGV >= 2){
    $newfile=$ARGV[2];
}else{
    $newfile=$fmof.".new";
}

print "New file is $newfile\n";

# GET CO2 coordinates
open(INc,"< $fco2");

$line_atom_number="    6    30    30    6\n";
$nOplus=12;
$nCplus=6;
$nOstart=24;	# 6(Zn) + 18(O)
$nCstart=60;	# 6(Zn) + 30(O) + 24(C)
$tag="OFF";
$natom_total=72;
$selective_dynamics="NO";
$natom=0;
@lattice_CO2=();
@lattice_MOF=();
@added_atom_kind=();
@added_atom_number=();
@added_coord=();
### read added atoms's coordinates
$i=0; $k=0; $l=0;
while($line=<INc>){
    next if($i<1);
    chomp($line);
    @field=split(/\s+/,$line);
    if($field[0] eq ""){ shift(@field);}
    if($i==1){ $added_scale=$field[0]; next;}
    if($i<=4){ 
	for($j=0;$j<=$#field;$j++){
#	    $lattice_CO2[$k]=$field[$j];
#	    $k++;
	    push(@lattice_CO2,$field[$j]);
	}   next;
    } 
    if($i==5){
	for($j=0;$j<=$#field;$j++){
            push(@added_atom_kind,$field[$j]);
	}   next;
    }
    if($i==6){
        for($j=0;$j<=$#field;$j++){
            push(@added_atom_number,$field[$j]);
        }   next;
    }
#    print @field[0]."\n";
    if($i==7){
	if($field[0] ne "Direct"){
	    print "ERROR: CO2 doesn't have direct coordinates"."\n";
	    exit;
	} next;
    }

#    print "$i= ", $tag,$natom,"\n";
    if($line ne ""){
	push(@added_coord,$line."\n");
	$l++;
    }else{next;}	
} continue {
    $i++;
}
close(INc); 

if($#added_atom_number != $#added_atom_kind){
    print "Error: different number of atom kinds in added_atom_kinds and added_atom_number\n";
}

if(sum(@added_atom_number) != $l){
    print "Error: different number of all atoms and number of line\n";
    exit;
}else{
    for($i=0,$j=0,$k=0;$i<=$#added_coord;$i++){
	print " Added: $added_atom_kind[$k]: $added_coord[$i]";
	if($j==$added_atom_number[$k]-1){
	    $j=0;
	    $k++;
	} $j++;
    }
}

$i=0;
$k_atom=-1;
$set_tag="OFF";
open(INm,"< $fmof");
open(OUT,"> $newfile");
@atoms_number=();
@atoms_number_orig=();
while($line=<INm>){
    if($i<5){ 
	print OUT $line;
	print $line;
	next; 
    }
    chomp($line);
    @field=split(/\s+/,$line);
    if($field[0] eq ""){ shift(@field);}

    if($i==5){
	# scan for all added atoms which should be somewhere in atom_list
	# index k for mother atoms and j for added atoms
	for($j=0;$j<=$#added_atom_kind;$j++){
	    $species_tag="F";
	    for($k=0;$k<=$#field;$k++){
		if($added_atom_kind[$j] eq $field[$k]){
		    $add_index[$j]=$k;
		    $species_tag="T";
		    next;
		}
	    }
	    if($species_tag eq "F"){
		print "Error: Added atoms is not listed in the mother poscar\n";
		exit;
	    }
	}
	print OUT $line."\n";
        print $line."\n";
	print "added atom indices: ", join("  ",@add_index)."\n";
	next;
    }
    if($i==6){ 
	for($k=0;$k<=$#field;$k++){
	    push(@atoms_number_orig,$field[$k]);
	    for($j=0;$j<=$#add_index;$j++){
		# if the number atom is added by add_index
		if($add_index[$j]==$k){
		    $field[$k]+=$added_atom_number[$j];
		}
	    }
	    push(@atoms_number,$field[$k]);
	}
	print OUT join("  ",@atoms_number)."\n";
	print     join("  ",@atoms_number)."\n";
	$set_tag="ON"; # start of $set_tag routine
	next;
    }

#   selective dynamics is not included at the moment
    next if($field[0] =~ /^[sS]/);
    if($field[0] =~ /^[dD]/) {
	print OUT $line."\n";
        print     $line."\n";
	next;
    }
    if($i_katom <= $katom_number-1){
#	print "i_katom $i_katom : katom_number $katom_number\n";
    	print OUT $line."\n";
	print     $line."\n";
    	$natom++;
	$i_katom++;
	# after adding i_katom i_katom == katom_number at the end
	if($i_katom == $katom_number){
	    # add coord line if there are added lines
	    if($add_tag eq "ON"){
	    	# the added coord should be added orderly
	    	for($al=0;$al<$added_number;$al++){
		    print OUT $added_coord[0];
		    print     "New: ",$added_coord[0];
		    $natom++;
		    shift(@added_coord);
	    	}
	    } 
	    $set_tag="ON";
    	}
    }
    
    if($natom == sum(@atoms_number)) {
	print "$natom lines are copied\n";
	last;
    }
} continue {
    $i++;
    # if ON, define katom_number [orig] and added_number for adding or not
    if($set_tag eq "ON"){
	# k_atom for species; i_katom is counting index
	$k_atom++;
	$katom_number=$atoms_number_orig[$k_atom];
	for($k=0;$k<=$#add_index;$k++){
	    if($k_atom == $add_index[$k]){
		$add_tag="ON";
		$added_number=$added_atom_number[$k];
		last;
	    }
	    $add_tag="OFF";
	}
	$i_katom=0;
 	$set_tag="OFF";
#	if($k_atom == 2){exit;}
    }
}
close(INm);
close(OUT); 
