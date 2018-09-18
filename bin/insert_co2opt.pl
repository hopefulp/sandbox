#!/usr/bin/perl
# written by Joonho Park
# read first CO2 coordinates as direct format and insert tham into mof POSCAR
# a.pl mof-file co2-file
if($#ARGV < 1){
    print "Error: input TWO POSCAR's ";
    exit(0);
}

$fin=$ARGV[0];
$fco2=$ARGV[1];
$new_suffix="add6co2";

if($#ARGV >= 2){
    $newfile=$ARGV[2];
}else{
    $newfile=$fin.".$new_suffix";
}

print "6CO2 is added in $newfile\n";

# GET CO2 coordinates
open(IN,"< $fco2");

$dynamics="Selective";
#$dynamics="Full Relaxation";

if($dynamics eq "Selective"){
    	$truefalse="  F  F  F\n";
}else{ 	$truefalse="  T  T  T\n";}

$i=0;
$iatom=0;
$ntatom_add=0;
@add_atom=();
@add_natom=();
### read three atoms's coordinates
LINE: while($line=<IN>){
    chomp($line);
    @field=split(/\s+/,$line);
    if($field[0] eq ""){ shift(@field);}

    next LINE if($i<5);
    if($i==5){ 
	for($j=0;$j<=$#field;$j++){
	    push(@add_atom,$field[$j]);
	} next LINE;
    }
    if($i==6){
        for($j=0;$j<=$#field;$j++){
	    $ntatom_add+=$field[$j];
            push(@add_natom,$field[$j]);
        } next LINE;
    }

    if($field[0] =~ /^[cC]/) {
        print "Error::The input file is expected to be Direct format\n";
	exit;
    }
    next LINE if($field[0] =~ /^[dD]/);
    
#    print "$iatom $ntatom_add\n";
    if($iatom < $ntatom_add){ 
	$coord[$iatom][0]=$field[0];
	$coord[$iatom][1]=$field[1]; 
	$coord[$iatom][2]=$field[2];
	$iatom++;	next LINE;
    }
} continue {
    $i++;
}
close(IN); 


#for($i=0;$i<$ntatom_add;$i++){
#    for($j=0;$j<3;$j++){
#    	print " $coord[$i][$j]";
#    } print "\n";
#}

$ikind=0;		# index of kind from 0 to 3 in the case of 4 atoms
$ikind_add=0;		# index for @add_atom() element
$iatom_kind=0;		# counting index for atoms in a certain kind
$iadd=0;		# counting index for O and C to be added
$i=0;
$test="OFF";		# tag for test-moment to check whether change the flow
$ntatoms=0;		# not used
@atom_kinds=();		# atom kinds for MOF, 
@natom_kind=();		# number of atoms in each kind
@ntatom_kind=();	# number of total atoms in each kind when add added-co2
open(IN,"< $fin");
open(OUT,"> $newfile");
LOOP: while($line=<IN>){
    chomp($line);
    @field=split(/\s+/,$line);
    if($field[0] eq ""){ shift(@field);}
    if($i<5) { 
	print $line,"\n";
	print OUT $line,"\n";
	next LOOP;
    }
    if($i==5){ 	# atom kinds
        if($field[0] !~ /[A-Z]/) {
            print "Error in $fin format: there are not atom symbols in 6th line\n";
            exit(0);
        }else{ for($j=0;$j<=$#field;$j++){push(@atom_kinds,$field[$j]);} }
	print join("  ",@atom_kinds)."\n";
        print OUT join("  ",@atom_kinds)."\n";
	next LOOP;
    }
    if($i==6){ 	# number of atoms in each species
        for($j=0;$j<=$#field;$j++){
	    $natom_kind[$j]=$field[$j];		# of atoms in each kind of bare MOF
	   # print "$j : $natom_kind[$j]\n";
	    $tag="OFF";
	    for($k=0;$k<=$#add_atom;$k++){
		if($add_atom[$k] eq $atom_kinds[$j]){
		    $tag="ON";
		    $ntatom_kind[$j]=$field[$j]+$add_natom[$k]; # total # of atoms in the kind
		    last;
		}
	    }
	    if($tag eq "OFF"){
		$ntatom_kind[$j]=$natom_kind[$j];
	    }
            $ntatoms+=$ntatom_kind[$j];
        }
        print join("  ",@ntatom_kind),"\n";
        print OUT join("  ",@ntatom_kind),"\n";
        next LOOP;
    }
# end of memory for atom species and numbers
    if($field[0] =~ /^[sS]/){
	print $line."\n";
	print OUT $line."\n";
	next LOOP;
    }
    if($field[0] =~ /^[cC]/) {
        print "Error::The input file is expected to be Direct format\n";
        exit;
    }
    if($field[0] =~ /^[dD]/){
	print $line."\n";
        print OUT $line."\n";
        next LOOP;
    }

    # write or skip test

    if($test eq "OFF" and $iatom_kind ==0){
	if($natom_kind[$ikind] == $ntatom_kind[$ikind]){
	    	$tag_add="OFF";
	}else{ 	$tag_add="ON";}
	$test="ON";
    }
    # start write input MOF coordinates by line
    if($test eq "ON"){
	if($iatom_kind < $natom_kind[$ikind]){
	    for($j=0;$j<3;$j++){ 
		print "  $field[$j]";
		print OUT "  $field[$j]";
	    }
	    print  $truefalse;
	    print OUT  $truefalse;
	    $iatom_kind++;
	    next LOOP if($iatom_kind != $natom_kind[$ikind]); 
	}
	# if there is another part save in the list, write here
	if($tag_add eq "ON"){
	    # count the previously added number in the @add_atom list
	    for($l=0,$sum_iadd=0;$l<$ikind_add;$l++){
		$sum_iadd+=$add_natom[$l];
	    }
	    for($k=0;$k<$add_natom[$ikind_add];$k++){
		$iadd=$sum_iadd+$k;
		for($j=0;$j<3;$j++){
                    print "  $coord[$iadd][$j]";
                    print OUT "  $coord[$iadd][$j]";
		}
           	print  "  T  T  T\n";
            	print OUT  "  T  T  T\n";
	    }
	    $ikind_add++;
	}
	$test="OFF";$iatom_kind=0; $ikind++; next LOOP;
    }
} continue {
$i++;
}
close(IN);
close(OUT); 
