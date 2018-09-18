#!/usr/bin/perl
# written by Joonho Park
# energy overlap between $atom1 $atom1l orbital and $atom2 $atom2l orbital
# energy is given by center and delta E or Emin and Emax

use lib '/qcfs/joonho/modules';
use VaspDos;

if($#ARGV<1){
    print "Usage :: $0 PROCAR [e=]-2.5[:E_max] \"atom list 3 9 etc\" \n";
    exit(0);
}



$switch_density="pdos";		# ldos for atom, pdos for angular momentum
$fin=$ARGV[0];

$str_energy=$ARGV[1];
if($str_energy =~ /=/ ){
    @e=split(/=/,$str_energy);
    if($e[1] =~ /:/){
        @f=split(/:/,$e[1]);
        $Emin=$f[0];
        $Emax=$f[1];
    }else{
	$E_pivot=$e[1];
    }
}else{ 
    $E_pivot=$str_energy;
}

if (! defined $Emin){   
    $del_energy=0.5;
    $Emin=$E_pivot-$del_energy;
    $Emax=$E_pivot+$del_energy;
}


print "Emin =  $Emin, Emax = $Emax\n";

@atoms=();
$str_atoms=$ARGV[2];
if($#ARGV == 2){
    @l_atoms=split(/\s+/,$str_atoms);
    if($l_atoms[0] eq ""){ shift @l_atoms;}
    for($i=0;$i<=$#l_atoms;$i++){
    	push(@atoms,$l_atoms[$i]);
    }
}else{
    print "Input atoms array in \"   \" ";
    exit(2);
    #@atoms=qw(1 3);		# for CO2 DeG-B4C
    #@atoms=qw(5 6);		# for DeG-OBC
    #@atoms=qw(1 25); 	# for MOF and CO2
}
$natoms=$#atoms+1;
print "atom list: ($natoms atoms) ", join("  ",@atoms), "\n";

##### 2. charge density for total density of atom (ldos), 0.02
#####                   for partial density of atom, less than 0.01
@density=qw(0.01 0.01 0.01 0.01 0.01 0.01 0.01 ); 	# 0.02 0.02 ; 0.002 0.002 ;	0.003 0.003 
#@density=qw(0.01 0.01 ); 	# 0.02 0.02 ; 0.002 0.002 ;	0.003 0.003 
$density_ldos=0.02;
$density_pdos=0.01;
$density_limit=$density_ldos;
if(! defined(@density)){
    @density=();
    for($i=0;$i<=$#atoms;$i++){
	push(@density, $density_limit);
    }
}

##### 1. if @atomsl is commented, it calculates total density (ldos) of the atom
#@atomsl=qw(2 1);
@pdos=();
@ldos=();
@m_orbital=(['s'],['py','pz','px'],['dxy','dyz','dz2','dxz','dx2']);

open(IN,"<$fin");

$i=0;
$iatom=0;
$sum=0;
$band=0;
$ncheck_energy=0;
$skip_tag="OFF";	# do not skip now
while($line=<IN>){
    next if($i<1);
    chomp($line);
    ### get important parameters from 2nd line
    if($i==1){
	@line=split(/\W+/,$line);
	if($line[0] eq ""){shift(@line);}
	$tot_kpoints=$line[3];
	$tot_bands=$line[6];
	$tot_atoms=$line[9];
	$max_skip=$tot_atoms+4;
#	print $tot_kpoints, $tot_bands, $tot_atoms, "\n";
   } 

    next if($skip_tag eq "ON");
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}
	

    if($line[0] eq "k-point"){
	$kpoint=$line[1];
    }
    if($line[0] eq "band"){
	$band=$line[1];
	if($line[3] eq "energy"){
	    $energy=$line[4];
	    if($Emin <= $energy and $energy <= $Emax){
		$skip_tag="OFF";	# execute the following sentences
	    }else{
		$skip_tag="ON";		# go to next until meet $energy
		$nskip=0;
	    }
	    $ncheck_energy++;
	    next;
	}else{ print "Error in getting energy\n"; exit(0);
	}
    }

    ############### following #########################
    ### only for the selected energy
    ### for the specific atoms, check the density populated at the same energy

    for($j=0;$j<=$#atoms;$j++){
    	if($line[0] == $atoms[$j]){
#   	    if(defined(@atomsl) and $ncheck_energy==1){	($nc,$nd)=&VaspDos::lm_quantum($atomsl[$j]);}
   	    if(defined(@atomsl)){	@column_list=&VaspDos::lm_quantum($atomsl[$j],$m,1); }
	    else{push(@column_list,9);} # this include only orbitals by "shift(@list); in sub sum()"
	    $nc=$column_list[0]; $nd=$column_list[$#column_list];
#	    print "index: ",$nc,"\t",$nd,"\n";
#	    print $j,$atomsl[$j],"\t",$nc,$nd,"\n";
	    $pdos[$j]=&sum($nc,$nd,@line);
	    $ldos[$j]=$line[10];	# this include all the line column so ldos lies on 11-th column
    	}
    }
    if($line[0] == $tot_atoms){
#    	print  $energy,"\t", $kpoint,"\t",$band,"\t",join("  ",@pdos),"\t", join("  ",@ldos),"\n";
	if($pdos[0] >= $density[0] and $pdos[1] >= $density[1]){
    	    print  $energy,"\t", $kpoint,"\t",$band,"\t", $pdos[0],"\t", $pdos[1],"\n";
	}
    }
  
} continue {
    $i++;
    ### whenever $tag turns ON $iatom is increased
    if($tag eq "ON"){
	$iatom++;
    }
    ### increase the number of skip
    if($skip_tag eq "ON"){
	$nskip++;
	### skip is over : turn off $skip_tag
	if($nskip >= $max_skip){
	    $skip_tag="OFF";
	}
    }
}
close(IN);
# summation in the dos line of PROCAR
sub sum{
    my($i,$f,@list)=@_;
    my($j,$sum);

    shift(@list);	# remove atom index
    $sum=0;
    for($j=$i;$j<=$f;$j++){
	$sum+=$list[$j];
#	print $list[$j],"\t";
    }
#    print "\n";
    return $sum;
}


