#!/usr/bin/perl
# written by Joonho Park
# read  DOSCAR or dir/DOSCAR if there is directory name
# read l and m quantum number
# for spin up or down s=[1|0] 1 for up and 0 for down
# for atom a:1:2 or a:25 or a=1:25

$suffix="";
$suffixl="";
$suffixlm="";
use lib '/home/joonho/modules';
use GetHash;
use VaspDos;

########## argument analysis 
@a_list=();
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
    	    }else{ push(@a_list,$num);}
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

print "quantum number: \t";
if(defined($l)){ print "l = $l\t";}
else { print "NA\t";}
if(defined($m)){ print "m = $m\n";}
else { print "NA\n";}

print "atom list: ",join(' ',@a_list),"\n";

######## ordering atoms and output file naming

@orbital=('s','p','d','f');
#@magnetic_orbital=(['y','z','x'],['xy','yz','z2','xz','x2-y2']); VASP-DOSCAR ordering
@magnetic_orbital=(['0'],['x','y','z'],['xy','yz','xz','x2-y2','z2']);
if(defined($l)){
    $suffix="_$orbital[$l]";
    $suffixl=$suffix;
    if(defined($m)){
	$suffix.=$magnetic_orbital[$l][$m];    
	$suffixlm=$suffix;
    }
}
if( @a_list ne "" ){
    if($#a_list == 0){
	$suffixa="_a$a_list[0]";
    }else{
    	$suffixa="_a$a_list[0]-$a_list[$#a_list]";
    }
    $suffix.=$suffixa;
    $suffixl.=$suffixa;
    $suffixlm.=$suffixa;
}
if(defined($spin)){
    $suffix.="_$spin";
    $suffixa.="_$spin";
    $suffixl.="_$spin";
    $suffixlm.="_$spin";
}

$fouta="Ldos$suffixa.dat";

### TDOS is overall
### LDOS, PDOS is spin up and down divided
print "Input: ",$fin."\t";
$fout="Tdos.dat";
$foutla="Pdos$suffixl.dat";
$foutlma="Pdos$suffixlm.dat";
print "OUT: Tdos $fout\t";
if(defined($l)){
    print "; Pdos $foutla\t";
    open(OUTla,"> $foutla");
}
if(defined($m)){
    print "; Pdos lm  $foutlma\t";
    open(OUTlma,"> $foutlma");
}

if(@a_list ne ""){
    print "; Ldos $fouta\t";
    open(OUTa,"> $fouta");
}

print "\n";

open(IN,"< $fin");
open(OUT,"> $fout");


$i=0;
$iatom=0;	# note first atom 0 is for total DOS not atom
$iene=0;
$Nene=1000;	# a large number to skip $Nene == $iene == 0
@energy=();
$tag_spin="NO";
$delimit_atom="";
$tag_delimit="OFF";
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
#    if($tag_delimit eq "OFF" && $iene==0 && $line[0] != -20.000){
    if($tag_spin eq "YES" and $tag_delimit eq "OFF" and $line eq $delimit_atom){
   	$tag_delimit="ON";	
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
    if($tag_delimit eq "ON"){
	$tag_delimit="OFF";
    }
#    print $iene,"\t", $Nene,"\n";
    if( $iene == $Nene){
#	print $i."\n";
	$iene=0;
	$iatom++;
    }
}

# iatom start from 0, TDOS then, 54 atoms' PDOS, and last added +1 so it should be natom+1
if($natom != $iatom-1){
    print "ERROR: in counting atoms $natom != $iatom-1\n";
    exit(0);
}
############################ END: load DOSCAR 
############################ shuffle DOS
# other options
$ene_Fermi = 0;		# do not shift w.r.t Fermi level
#if(!defined($spin)){ $nspin="YES";}  # For DOSCAR, for PROCAR use "NO"
$nspin=2; 	#For DOSCAR, for PROCAR use 1
#($nc,$nd,$mc,$md)=&VaspDos::lm_quantum($l,$m,$nspin,$spin);
(@nl_col)=&VaspDos::lm_quantum($l,$m,$nspin,$spin);
print "column list : ", join(" ",@nl_col),"\n";

for($i=0;$i<$Nene;$i++){
#    if($energy[$i] == 3.013) { print $tdos[0][$i][0],  $tdos[0][$i][1],  $tdos[0][$i][2],  $tdos[0][$i][3],"\n"; exit; }
    printf "TDOS: Emin = %10.4f  Emax= %10.4f\n",$Emin, $Emax if($i==0);
    if($Lspin eq "F"){ # this is not complete
	# for TDOS, $j==$
    	print $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0],"\n";
    	print OUT $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0],"\n";
    }else{
	# SPIN
	# TDOS
#	print $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0]+$tdos[0][$i][1],"\n";
	print OUT $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0]+$tdos[0][$i][1],"\n";
	# for LDOS
	$pdos=0;
	$pdosm=0;
	$ldos=0;
	# LDOS if atoms are selected, LDOS is spin adapted if $spin is defined
 	if(@a_list ne ""){
	    # sum all the selected atoms	j atom index
	    for($j=0;$j<=$#a_list;$j++){	
		# sum all the orbitals; 	k dos colume index scans all the column
		for($k=0;$k<$ncol;$k++){
		    if(defined($spin) and $spin =~ /[dD]/ and $k%2==0){
			if($i==0) {print "$k is skipped for Down spin\n";}
			next;
		    }	# upspin is skipped
		    if(defined($spin) and $spin =~ /[uU]/ and $k%1==0){
			if($i==0) {print "$k is skipped for Down spin\n";}
			next;
		    }
		    $ldos += $tdos[$a_list[$j]][$i][$k];
		    # PDOS if orbitals are selected, PDOS is spin adapted
	    	    if(defined($l) && $nc <= $k && $k<=$nd){
		    	$pdos += $tdos[$a_list[$j]][$i][$k];
			if($mc <= $k && $k <= $md){
			    $pdosm += $tdos[$a_list[$j]][$i][$k];
			}
		    }
		}
	    }
		
	 # if atoms are not selected
	}else{
	    # PDOS for all atoms
	    if(defined($l)){
	    # sum all the atoms
	    	for($j=0;$j<=$#a_list;$j++){
	            for($k=$nc;$k<=$nd;$k++){
	        	$pdos += $tdos[$a_list[$j]][$i][$k];
	        	}
	        }
	    }
    	}
#	print $energy[$i],"\t", $pdos,"\n";
#	print $energy[$i]-$ene_Fermi,"\t", $pdos,"\n";
	if(defined($m)) {print OUTlma $energy[$i]-$ene_Fermi,"\t", $pdos,"\n";}
	if(defined($l)) {print OUTla  $energy[$i]-$ene_Fermi,"\t", $pdos,"\n";}
	print OUTa   $energy[$i]-$ene_Fermi,"\t", $ldos,"\n";
    }
}

close(IN); close(OUT); 
if(defined($l)){
    close(OUTla);
}
if(defined($m)){
    close(OUTlma);
}
if(0<=$#a_list){
    close(OUTa);
}

