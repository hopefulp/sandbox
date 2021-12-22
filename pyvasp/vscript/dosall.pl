#!/usr/bin/perl
# written by Joonho Park
# ? read  DOSCAR or dir/DOSCAR if there is directory name
# cd to inside of directory and run ./dosall.pl
### usage: ./dosall.pl a=1 l=1 m=0 s=1 or -1
# read l and m quantum number
# for spin up or down s=[1|0] 1 for up and 0 for down
# for atom  a=25 or a=1:25

$suffix="";
use lib '/home/joonho/sandbox/mod_perl';
use GetHash;
use VaspDos;

$Fermi_shift='F';	# T for shift, F for no-shift
$Fermi_check=0;		# write Fermi.dat at value EF with value 1.2 * max(Tdos)
$Max_Tdos=0;		# control of y-axis of Dos.dat

#print "Now run in the directory\n";

########## argument analysis 
@a_list=();

if($#ARGV < 0){
    print "Usage:: $0 [DOSCAR] for Tdos and Fermi.dat \n";
    print "Usage:: $0 a=na[:nb] l=l m=m s=[1 for up|-1 for down]\n";
    print "        below a are options\n";
    exit(1);
}

$fin="DOSCAR";
if($ARGV[0] =~ /^d/){
    @dir=split(/=/, $ARGV[0]);
    $fin="$dir[1]/DOSCAR";
    shift @ARGV;
    print "Note:: The outfiles will be located at present directory\n";
}
#exit(1);

for($i=0;$i <= $#ARGV; $i++){
    $atom_list="";
    if($ARGV[$i] =~ /=/){
	($quantum,$num)=&GetHash::getvalue($ARGV[$i],'=');
	# usage: l=2 m=2 or just l=1 then m is summed
        if($quantum eq "l"){
            $l=$num;
        }elsif($quantum eq "m"){
            $m=$num;
        }elsif($quantum eq "s"){
            $spin=$num;
        }elsif($quantum eq "a"){
	    $atom_list=$num;
    	    if($atom_list =~ /:/){
	    	    ($a_start,$a_end)=&GetHash::get_2values($atom_list,':');
	    	    for($j=$a_start;$j<=$a_end;$j++){
	    	        push(@a_list,$j);
	    	    }
    	    }elsif($atom_list =~ /,/){
                @a_list=split(', ', $atom_list)
            }else{ push(@a_list,$num);}
	}elsif($quantum eq "dir"){
	    $fin="$num/DOSCAR";
	}
    # if there is not non-digital (only digital), it's l or m
    }elsif($ARGV[$i] !~ /\D/){
	if(!defined($l)){ $l=$ARGV[$i];}
	elsif(!defined($m)){ $m=$ARGV[$i];}
    # if there is not non-word(a-Z,0-9,_) character, it's input file
    }elsif($ARGV[$i] !~ /\W/){
	if($ARGV[$i] eq "DOSCAR") {$fin = "DOSCAR";}
	else {$fin = "$ARGV[$i]/DOSCAR";}
    }
}

if(!defined($fin)) {$fin="DOSCAR";}
print "input = $fin\n";
print "quantum number: \t";
if(defined($l)){ print "l = $l\t";}
else { print "NA\t";}
if(defined($m)){ print "m = $m\t";}
else { print "NA\t";}
if(defined($spin)){ print "s = $spin\n";}
else { print "NA\n";}

print "atom list: ",join(' ',@a_list),"\n";

######## ordering atoms and output file naming

@orbital=('s','p','d','f');
#@magnetic_orbital=(['y','z','x'],['xy','yz','z2','xz','x2-y2']); VASP-DOSCAR ordering
@magnetic_orbital=(['0'],['x','y','z'],['xy','yz','xz','x2-y2','z2']);
if(defined($l)){			# if $l is defined, make Pdos
    $suffix="_$orbital[$l]";
    if(defined($m)){
	$suffix.=$magnetic_orbital[$l][$m];    
    }
    $fdos="Pdos";
}

if( 0 <= $#a_list ){
#    print "atom_list: ",$#a_list, join("  ",@a_list),"\n";
    if($#a_list == 0){
	$suffix.="_a$a_list[0]";
    }else{
    	$suffix.="_a$a_list[0]-$a_list[$#a_list]";
    }
    if(!defined($fdos)){ $fdos="Ldos";}	# if $a is defined but $l is not defined make Ldos
}

if(defined($spin)){
    if($spin==1){ $spin_l="up";}
    if($spin==-1){ $spin_l="down";}
    $suffix.="_$spin_l";
}

if(!defined($fdos)){ $fdos="Tdos";}

$fout="$fdos$suffix.dat";
if($Fermi_shift eq "T"){
    $fout="$fdos$suffix"."F0.dat";
}

### TDOS is overall
### LDOS, PDOS is spin up and down divided
print "Input: $fin\t";
print "OUT: $fout\n";

open(IN,"< $fin");
open(OUT,"> $fout");
if($Fermi_shift eq "F" and $#a_list < 0){
    open(EFOUT,"> Fermi.dat");
}

$i=0;
$iatom=0;	# note first atom 0 is for total DOS not atom
$iene=0;
$Nene=1000;	# a large number to skip $Nene == $iene == 0
@energy=();
$tag_spin="NO";
$delimit_atom="";
$Tag_Atom_delimit="OFF";
###################################################  read three atoms's coordinates
LINE: while($line=<IN>){
#    if($iene==0) {print $line;}
#    if($line eq $delimit_atom){ next LINE;}
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($i==0){ $natom=$line[0]; next LINE;}
    if($i==1){ next LINE;} # add something here if you know DOSCAR
    if($i==2){ next LINE;}
    if($i==3){ $threeletter=$line[0]; next LINE;}
    if($i==4){ $model=$line[0]; next LINE;}
    if($i==5){
	$delimit_atom=$line;
	$Emax=$line[0];
	$Emin=$line[1];
	$Nene=$line[2];
	$ene_Fermi=$line[3];
	$dnknow=$line[4];
	next LINE;
    }
    # this part locates after $i==6  for $iene=0 and start in an atom block 
    if($line eq $delimit_atom){
   	$Tag_Atom_delimit="ON";	
	next LINE;
    }


### DOS LINE: 1. energy 2. tdos_@_E    3. tdos_int_E 						for non-spin
### DOS LINE: 1. energy 2. tdos_@_E_up 3. tdos_@_E_down 4. tdos_int_E_up 5. tdos_int_E_down 	for spin
    # define spin 
    if($tag_spin eq "NO"){
        if($#line == 2){
            $Lspin='F';
        }elsif($#line == 4){
            $Lspin='T';
        }else{ print "Spin is not determined so die\n"; exit(1);} 
	$tag_spin="YES";
    }
    # save all the components in dos line
#    print "iatom = $iatom\n";
    # $iatom == 0 is for total DOS; atoms start from 1 to # of atoms
    if($iatom == 0){ # for total
	push(@energy,$line[0]);
#    	print $iatom,"\t",$line[0],$energy[$iene],"\n";
    }
    # save all the dos components: first TDOS then atoms
    for($j=1;$j<=$#line;$j++){
        $tdos[$iatom][$iene][$j-1]=$line[$j];
	if($iatom == 1){$ncol=$#line;} # number of dos component
    }
    $iene++; #print "iene = $iene\n";

} continue {
    $i++;
    # when you read delimiter between PDOS block, skip to next line withough anithing changed, $iene, $iatom
    # So the number of atom is same as the number of delimit
    if($Tag_Atom_delimit eq "ON"){
	$Tag_Atom_delimit="OFF";
	if($iatom==1){ $Max_iene=$iene;}
	$iene=0;
	$iatom++;
    }
}

# iatom start from 0, TDOS then, 54 atoms' PDOS, and last added +1 so it should be natom+1
if($natom != $iatom){
    print "ERROR: in counting atoms $natom != $iatom\n";
    exit(0);
}
############################ END: load DOSCAR 
############################ shuffle DOS
#if(!defined($spin)){ $nspin="YES";}  # For DOSCAR, for PROCAR use "NO"
#($nc,$nd,$mc,$md)=&VaspDos::lm_quantum($l,$m,$nspin,$spin);
#### other options

$Efermi=$ene_Fermi;

if($Fermi_shift eq "F"){
    $ene_Fermi = 0;		# do not shift w.r.t Fermi level
    $EF_test=$Efermi;
    print "Energy is not shifted\n";
}else{
    print "Energy is shifted to Fermi Level\n";
    $EF_test=0;
}

if($Lspin eq "T"){
    $nspin=2;
}elsif($Lspin eq "F"){
    $nspin=1;
}else{
    print "Magnetism error\n";
    exit(5);
}

#For DOSCAR, for PROCAR use 1 , note: this is not for spinless
#column list start from 0 to 17 in the case of s(2), p(6), d(10) for spin
# For non-magnetic 
(@nl_col)=&VaspDos::lm_quantum($l,$m,$nspin,$spin);
print "column list : ", join(" ",@nl_col),"\n";
for($i=0;$i<$Nene;$i++){
    if($Max_Tdos < $tdos[0][$i][0]){
	$Max_Tdos = $tdos[0][$i][0];
    }
}
print "Max_Tdos = $Max_Tdos\n";
for($i=0;$i<$Nene;$i++){
#    if($energy[$i] == 3.013) { print $tdos[0][$i][0],  $tdos[0][$i][1],  $tdos[0][$i][2],  $tdos[0][$i][3],"\n"; exit; }
    printf "TDOS: Emin = %10.4f  Emax= %10.4f Fermi level= %10.4f\n",$Emin, $Emax, $Efermi if($i==0);
    ### This is for non-magnetic DOSCAR
    if($Lspin eq "F"){ 
    	#print $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0],"\n";
    	### tdos[0] is for TDOS and others for atoms
    	### tdos[][][0] is for DOS at the energy ; tdos[][][1] is for accumulated DOS
        ### if there is no atom list, write Tdos
        if($#a_list < 0){
           print OUT $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0],"\n";
           if($Fermi_shift eq "F" and $#a_list < 0){
               if($Fermi_check==0 and $energy[$i] < $EF_test){
                   print EFOUT $energy[$i]-$ene_Fermi,"\t 0\n";
               }elsif($Fermi_check==0 and $energy[$i] >= $EF_test){
                   print EFOUT "$EF_test\t", $Max_Tdos*2, "\n";
                   $Fermi_check=1;
               }else{
                   print EFOUT $energy[$i]-$ene_Fermi,"\t 0\n";
               }
           }
        ### For atom list, get LDOS or PDOS depending on l quantum number
        }else{
	        $pdos=0;
            $pdosm=0;
            $ldos=0;
            ### as for an energy, loop atom lists [whethe not sorted], add l columns
            for($j=0;$j<=$#a_list;$j++){
                for($ic=0;$ic<=$#nl_col;$ic++){
                    $k=$nl_col[$ic];
                    $ldos += $tdos[$a_list[$j]][$i][$k];
                    #print "ldos $ldos\n";
                    }
                }
            print OUT   $energy[$i]-$ene_Fermi,"\t", $ldos,"\n";
	    }

    ### This is for magnetic DOSCAR
    }else{
#	print $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0]+$tdos[0][$i][1],"\n";
#	### if there is no atom list, write Tdos
	if($#a_list < 0){
	    print OUT $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0]+$tdos[0][$i][1],"\n";
	    if($Fermi_shift eq "F" and $#a_list < 0){
	        if($Fermi_check==0 and $energy[$i] < $EF_test){
    	            print EFOUT $energy[$i]-$ene_Fermi,"\t 0\n";
	        }elsif($Fermi_check==0 and $energy[$i] >= $EF_test){
	            print EFOUT "$EF_test\t", $Max_Tdos*2, "\n";
	            $Fermi_check=1;
	        }else{
	            print EFOUT $energy[$i]-$ene_Fermi,"\t 0\n";
	        }
	    }
	### For atom list, get LDOS or PDOS depending on l quantum number
	}else{
	    # for LDOS
	    $pdos=0;
	    $pdosm=0;
	    $ldos=0;
	    # LDOS if atoms are selected, LDOS is spin adapted if $spin is defined
 	
	    # sum all the selected atoms	j atom index
	    for($j=0;$j<=$#a_list;$j++){	
		#print "sum $a_list[$j] -th atom\n";
		#print join(" ",@nl_col),"\n";
		# sum all the orbitals; 	k dos colume index scans all the column
		for($ic=0;$ic<=$#nl_col;$ic++){
		    $k=$nl_col[$ic];
		    $ldos += $tdos[$a_list[$j]][$i][$k];
		    #print "ldos $ldos\n";
		}
	    }
	    print OUT   $energy[$i]-$ene_Fermi,"\t", $ldos,"\n";
		
	 # if atoms are not selected
#	}else{
#	    # PDOS for all atoms
#	    if(defined($l)){
#	    # sum all the atoms
#	    	for($j=0;$j<=$#a_list;$j++){
#	            for($k=$nc;$k<=$nd;$k++){
#	        	$pdos += $tdos[$a_list[$j]][$i][$k];
#	        	}
#	        }
#	    }
    	}
#	print $energy[$i],"\t", $pdos,"\n";
#	print $energy[$i]-$ene_Fermi,"\t", $pdos,"\n";
    }
}

close(IN); close(OUT); 

