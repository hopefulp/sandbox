#!/usr/bin/perl
# written by Joonho Park
# read mol file in the same lattice 
# only direct coordinate
# a.pl bare_mof mol_in  "C 2 H 2"

if($#ARGV < 2 or $ARGV[0] =~ /Usage/i){
    print "Error: input TWO POSCAR's and optional [outfilename] MOL TYPE \n";
    print "Usage: $0 bare_mof mol_in_lattice \"atom number\"\n";
    exit(1);
}#elsif($#ARGV > 2){
#    print "Error: Usage: Mol Type should be separate with blank quoted by \" and \" \n";
#    exit(0);
#}

$fmof_emp=$ARGV[0];
$fmol_lattice=$ARGV[1];
$new_suffix="addmol";

$newfile=$fmof_emp.".$new_suffix";
$add_mol=$ARGV[2];

print "OUTPUT: $newfile\n";

if($add_mol=~/ /){
    @add_mol=split(/ +/,$add_mol);
}else{
    print "Error: Usage: Mol Type should be separate with blank quoted by \" and \" \n";
    exit(2);
}

while(@add_mol){$add_kind[$i]=shift(@add_mol); $add_num[$i]=shift(@add_mol); $i++;}

for($i=0;$i<=$#add_kind;$i++){
    print "Atom: $add_kind[$i] Number: $add_num[$i] \n";
}

print "MOL will be added in $newfile\n";

# GET MOL coordinates
open(IN,"< $fmol_lattice");

$dynamics="Selective";
#$dynamics="Full Relaxation";

if($dynamics eq "Selective"){
    	$truefalse="  F  F  F\n";
}else{ 	$truefalse="  T  T  T\n";}

$i=0;
$iatom=0;
$nt_add=0;
@add_atom_kind=();
@add_atom_coord=();
@add_atom_num=();

### read MOL's coordinates
print "mol + mof contains all kinds of atoms\n";
LINE: while($line=<IN>){
    next LINE if($i<5);

    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){ shift(@line);}

    if($i==5){
	$string=join("",@line);
	die "There is numeric in atom line $!\n" if($string =~ /[^a-zA-Z]/); # /\d/
	for($j=0;$j<=$#line;$j++){
	    push(@add_atom_kind,$line[$j]);
	} next LINE;
    }
    if($i==6){
	$string=join("",@line);
	die "There is non-numeric in number line $!\n" if($string =~ /\D/);
        for($j=0;$j<=$#line;$j++){
	    $nt_add+=$line[$j];
            push(@add_atom_num,$line[$j]);
        }
    }
    if($i==7){
        if($line[0] =~ /^[cC]/) {
            print "Error::The input file is supposed to be Direct coordinates\n"; exit 11;
        }elsif($line[0] =~ /^[dD]/){
            $new_line[$i]=$line."\n";
            print "Direct coordinate is OK.\n"; next LINE;
        }else{ print "format error in cC|dD with $line[0]\n"; exit 12;}
    }

    # nt_add is calculated before
    if($i>=8){
        $add_atom_coord[$iatom]=$line."\n";
        $iatom++;
        if($iatom<$nt_add) {next LINE;}
        else {last LINE;}
    }

} continue {
    $i++;
}
close(IN); 


$i=0;
@bare_mof_coord=();		# atom kinds for MOF, 
@bare_mof_kind=();		# number of atoms in each kind
@bare_mof_num=();
$bare_nt_atom=0;
$acc_natom=0;		# accumulated total atom index in original MOF
$acc_l=0;
$iatom=0;		# index for total atom of original MOF
$nt_atom=0;		# index for atom species of original MOF
$insert_tag='F';
open(IN,"< $fmof_emp");
LOOP: while($line=<IN>){
    if($i<5) {  next LOOP; }
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){ shift(@line);}
    if($i==5){
        $string=join("",@line);
        die "There is numeric in atom line $!\n" if($string =~ /[^a-zA-Z]/);
        for($j=0;$j<=$#line;$j++){push(@bare_mof_kind,$line[$j]);}
        #$new_line[$i]=$new_line_kind."\n";
        next LOOP;
    }
    if($i==6){
        for($j=0;$j<=$#line;$j++){
            $bare_nt_atom+=$line[$j];
            $bare_mof_num[$j]=$line[$j];   # number of atoms of each kind
        }
        #$new_line[$i]=$new_line_num."\n";
        next LOOP;
    }
    next LOOP if($line[0] =~ /^[sS]/);  # for selected dynamics

    if($i==7){
        if($line[0] =~ /^[cC]/) {
            print "Error::The input file is supposed to be Direct coordinates\n"; exit 11;
        }elsif($line[0] =~ /^[dD]/){
            $new_line[$i]=$line."\n";
            print "Direct coordinate is OK.\n"; next LOOP;
        }else{ print "format error in cC|dD with $line[0]\n"; exit 13;}
    }

    # ntotal is calculated before
    if($i>=8){
        $bare_mof_coord[$iatom]=$line."\n";
        $iatom++;
        if($iatom<$bare_nt_atom) {next LOOP;}
        else {last LOOP;}
    }


} continue {
    if($i==0){  $new_line[$i]=join("",@add_atom_kind); }
    if($i<5) {  $new_line[$i]=$line;}
    $i++;
}
close(IN);

#### atom kind arrangement
%periodic_table=( H => 1, C => 6, N => 7, O => 8,
		Mg => 12, Ti => 22, V => 23, Mn => 25, Fe => 26, 
		Co => 27, Ni => 28, Cu => 29, Zn => 30); 
@new_kind=();
@new_num=();
@new_coord=();
@pivot=();
### reference to MOF, check and add mol
for($i=0,$j=0,$iatom=0;$i<=$#bare_mof_kind;$i++){
    print $bare_mof_kind[$i],"\t",$add_atom_kind[$j],"\n";
    if($bare_mof_kind[$i] eq $add_atom_kind[$j]){
	push(@new_kind,$bare_mof_kind[$i]);
	$new_nk_atom=$bare_mof_num[$i]+$add_atom_num[$j];
	push(@new_num,$new_nk_atom);
	$j++;
	push(@pivot,2);
    }else{
	# only MOF
	print $periodic_table{$bare_mof_kind[$i]},"\t",$periodic_table{$add_atom_kind[$j]},"\n";
	if($periodic_table{$bare_mof_kind[$i]} > $periodic_table{$add_atom_kind[$j]}){
	    push(@new_kind,$bare_mof_kind[$i]);
	    push(@new_num,$bare_mof_num[$i]);
	    push(@pivot,1);
	# only mol
	}elsif($periodic_table{$bare_mof_kind[$i]} < $periodic_table{$add_atom_kind[$j]}){
	    push(@new_kind,$add_atom_kind[$j]);
	    push(@new_num,$add_atom_num[$j]);
	    $j++;
	    push(@pivot,-1);
	    redo;
	}else{	print "Error in merging atoms\n"; exit 10; }
    }
}

### if mol file has extra atom kinds
if($j<=$#add_atom_kind){
    while($j<=$#add_atom_kind){
	push(@pivot,-1);
	$j++;
    }
}

print join(" ",@pivot),"\n";
open(OUT,"> $newfile");
for($i=0; $i < 8; $i++){
    if($i==5){ 
	print join("  ",@new_kind)."\n";
	print OUT join("  ",@new_kind)."\n";
    }elsif($i==6){
        print join("  ",@new_num)."\n";
        print OUT join("  ",@new_num)."\n";
    }else{
    	print $new_line[$i];
    	print OUT $new_line[$i];
    }
}
### i for mof, j mol, k all, l=atom
for($i=0,$j=0,$k=0; $k<=$#new_kind; $k++){
    if($pivot[$k]==2){
	for($l=0;$l<$bare_mof_num[$i];$l++){
	    $string = shift(@bare_mof_coord);
	    print $string; print OUT $string;} $i++;
	for($l=0;$l<$add_atom_num[$j];$l++){ 
            $string = shift(@add_atom_coord);  
            print $string; print OUT $string;} $j++;
    }elsif($pivot[$k]==1){
	for($l=0;$l<$bare_mof_num[$i];$l++){ 
	    $string = shift(@bare_mof_coord);
	    print $string; print OUT $string; 
	} $i++;
    }elsif($pivot[$k]==-1){
	for($l=0;$l<$add_atom_num[$j];$l++){ 
	    $string = shift(@add_atom_coord);
	    print $string; print OUT $string; 
	} $j++;  
    }else {print "Error in printing\n";}
}
close(OUT); 

