#!/usr/bin/perl
# written by Joonho Park
# read first original MOF then
# MOF-CO2 files to get CO2 coordinates as direct format and insert them into original MOF
# a.pl mof1 mof2_CO2 [outf=filename] "C 2 H 2"

if($#ARGV < 2){
    print "Error: input TWO POSCAR's and optional [outfilename] MOL TYPE \n";
    exit(0);
}#elsif($#ARGV > 2){
#    print "Error: Usage: Mol Type should be separate with blank quoted by \" and \" \n";
#    exit(0);
#}

$fmof1=shift(@ARGV);
$fmof2=shift(@ARGV);
$new_suffix="addmol";

$tmp=shift(@ARGV);
if($tmp=~/=/){
    @fname=split(/=/,$tmp);
    $newfile=$fname[1];
    $mol=shift(@ARGV);
}else{
    $newfile=$fmof1.".$new_suffix";
    $mol=$tmp;
}

print "OUTPUT: $newfile\n";

if($mol=~/ /){
    @mol=split(/ +/,$mol);
}else{
    print "Error: Usage: Mol Type should be separate with blank quoted by \" and \" \n";
    exit(0);
}

$i=0;
while(@mol){$add_kind[$i]=shift(@mol); $add_num[$i]=shift(@mol); $i++;}

for($i=0;$i<=$#add_kind;$i++){
    print "Atom: $add_kind[$i] Number: $add_num[$i] \n";
}

print "MOL will be added in $newfile\n";

# GET MOL coordinates
open(IN,"< $fmof2");

$dynamics="Selective";
#$dynamics="Full Relaxation";

if($dynamics eq "Selective"){
    	$truefalse="  F  F  F\n";
}else{ 	$truefalse="  T  T  T\n";}

$i=0;
$l=0;
$iatom=0;
$iadded=0;
$ntatom=0;
@mof_katom=();
@mof_knatom=();
@extr_index=();

### read MOL's coordinates
LINE: while($line=<IN>){
    next LINE if($i<5);

    chomp($line);
    @field=split(/\s+/,$line);
    if($field[0] eq ""){ shift(@field);}

    if($i==5){ 
	for($j=0;$j<=$#field;$j++){
	    push(@mof_katom,$field[$j]);
	} next LINE;
    }
    if($i==6){
        for($j=0;$j<=$#field;$j++){
	    $ntatom+=$field[$j];
            push(@mof_knatom,$field[$j]);
        } 
### Extract MOL
	$acc_natom=0;
	for($j=0;$j<=$#mof_katom;$j++){
	    $acc_natom+=$mof_knatom[$j];
	    for($k=0;$k<=$#add_kind;$k++){
		if($mof_katom[$j] eq $add_kind[$k]){
		    $i_index=$acc_natom-$add_num[$k]+1;
		    $f_index=$acc_natom;
#		    print "i= $i_index, f= $f_index\n";
		    @p_index=($i_index..$f_index);
		    push(@ext_index,@p_index);
		    last;
		}
	    }
	}
	next LINE;
    }
### skip word lines
    if($field[0] =~ /^[cC]/) {
        print "Error::The input file is supposed to be Direct format\n";
	exit;
    }
    next LINE if($field[0] =~ /^[dDsS]/);

### read atom coordinates
#    print "$iatom $ntatom\n";
    if($iatom==$ntatom) {last; } # print "$i-th line\n";
    if($field[0]!~/\d/){ print "Error: Here is not coordinate line in $i-th, $iatom iatom, $ntatom ntatom\n";}
    if($iatom+1 == $ext_index[$l]){ 
	$coord[$iadded][0]=$field[0];
	$coord[$iadded][1]=$field[1]; 
	$coord[$iadded][2]=$field[2];
	$l++;
	$iadded++;
    }
    $iatom++;
} continue {
    $i++;
}
close(IN); 

print join(" ",@ext_index),"\n";

for($i=0;$i<=$#ext_index;$i++){
    for($j=0;$j<3;$j++){
    	print " $coord[$i][$j]";
    } print "\n";
}

$ikind=0;		# index of kind from 0 to 3 in the case of 4 atoms
$i=0;
$ntatom2=0;		# not used
@atom_kinds=();		# atom kinds for MOF, 
@natom_kind=();		# number of atoms in each kind
@natom_kind_all=();
$acc_natom=0;		# accumulated total atom index in original MOF
$acc_l=0;
$iatom=0;		# index for total atom of original MOF
$ikatom=0;		# index for atom species of original MOF
$insert_tag='F';
open(IN,"< $fmof1");
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
            print "Error in $fmof1 format: there are not atom symbols in 6th line\n";
            exit(0);
        }else{ for($j=0;$j<=$#field;$j++){push(@atom_kinds,$field[$j]);} }
	print join("  ",@atom_kinds)."\n";
        print OUT join("  ",@atom_kinds)."\n";
	next LOOP;
    }
    if($i==6){ 	# number of atoms in each species
        for($j=0;$j<=$#field;$j++){
	    $natom_kind[$j]=$field[$j];		# of atoms in each kind of bare MOF
	    $natom_kind_all[$j]=$field[$j];
	    for($k=0;$k<=$#add_kind;$k++){
		if($add_kind[$k] eq $atom_kinds[$j]){
		    $natom_kind_all[$j]+=$add_num[$k]; # total # of atoms in the kind
		    last;
		}
	    }
            $ntatom2+=$natom_kind[$j];
        }
        print join("  ",@natom_kind_all),"\n";
        print OUT join("  ",@natom_kind_all),"\n";
	$acc_natom=$natom_kind[0]; # to write to the end of 1st species
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

    # Get the index for writing of original MOF

    if($iatom==$acc_natom){ # determine to insert atoms after writing atom coordinates of one species
	### matching test
	for($k=0;$k<=$#add_kind;$k++){
	    if($atom_kinds[$ikatom] eq $add_kind[$k]){
		$insert_tag='T';
		$insert_index=$k;
		last;
	    }
	}
	$ikatom++;	# go to the next species
	$acc_natom+=$natom_kind[$ikatom];
    }

    #### insert part for added atoms
    if($insert_tag eq 'T'){
	for($l=1;$l<=$add_num[$insert_index];$l++){
	    for($j=0;$j<3;$j++){
                print "  $coord[$acc_l+$l-1][$j]";
                print OUT "  $coord[$acc_l+$l-1][$j]";
            }
            print  "  T  T  T\n";
            print OUT  "  T  T  T\n";
	}
	$acc_l+=$add_num[$insert_index];
 	$insert_tag='F';
    }

    ### Last insert part is finished
    if($iatom==$ntatom2) {last;}

    #### go over to the next species
    for($j=0;$j<3;$j++){
    	print "  $field[$j]";
	print OUT "  $field[$j]";
    }
    print  $truefalse;
    print OUT  $truefalse;
    $iatom++;
} continue {
$i++;
}
close(IN);
close(OUT); 
