#!/usr/bin/perl
# written by Joonho Park
# energy overlap between $atom1 $atom1l orbital and $atom2 $atom2l orbital
# energy is given by center and delta E or Emin and Emax

use lib '/qcfs/joonho/modules';
use VaspDos;

if($#ARGV<1){
    print "It needs input file and energy for detection\n";
    exit(0);
}

$fin=$ARGV[0];

$pivot_energy=$ARGV[1];
$del_energy=0.5;
$Emin=$pivot_energy-$del_energy;
$Emax=$pivot_energy+$del_energy;
if($#ARGV == 2){
    $Emin=$pivot_energy;
    $Emax=$ARGV[2];
}

##### 1. if @atomsl is commented, it calculates total density (ldos) of the atom
#@atomsl=qw(2 2 1 1);

##### 2. charge density for total density of atom (ldos), 0.02
#####                   for partial density of atom, less than 0.01
#@density=qw(0.02 0.02);
#@density=qw(0 0 0 0 0 0 0.00 0.000);
@density=qw(0.003 0.003);

@atoms=qw(1 3);		# for CO2, DeG-B4C
#@atoms=qw(5 6); 	# for DeG-OBC
#@atoms=qw(1 2 3 4 5 6 25 26); 	# for MOF and CO2
@atom_kind=();
@atomsl=();
for($i=0;$i<=$#atoms;$i++){
    if($atoms[$i] <= 6){
	push(@atom_kind,"Me");
	push(@atomsl,2); 	# for d-orbital
    }else{
	push(@atom_kind,"O");
	push(@atomsl,1); 	# for p-orbital
    }
}

$chg_diff_OO=0.00;	# 0.01
$chg_diff_MeMe=0.000; 	# 0.001

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
$band_max=0;
$band_min=10000;
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
    ### analyze after reading a k and a band
    if($line[0] == $tot_atoms){
#    	print  $energy,"\t", $kpoint,"\t",$band,"\t",join("  ",@pdos),"\t", join("  ",@ldos),"\n";
	$diff_OCO=$pdos[6]-$pdos[7]; # how much decreased in O of CO2
	$ave=average($pdos[1],$pdos[2],$pdos[3],$pdos[4],$pdos[5]);
	$diff_me=$pdos[0]-$ave;
#	if($pdos[6] >= $density[7] and $pdos[6] >= $density[7] and $diff_OCO >= 0.002){
#  	if($band == 86){
#	if($pdos[6]>0 and $density[7]>0){
	if($pdos[6]>$density[6] or $pdos[7]>$density[7]){
	    if($band_max < $band){ $band_max=$band;}
	    if($band < $band_min){ $band_min=$band;}
#    	    printf "%10.4f  %2d  %5d  %8.4f  %8.4f  %8.4f    %8.4f  %8.4f  %8.4f\n", $energy,$kpoint,$band,$pdos[0],$ave,$diff_me,$pdos[6],$pdos[7],$diff_OCO;
    	    printf "%10.4f  %2d  %5d  %8.4f  %8.4f  %8.4f    %8.4f  %8.4f  %8.4f\n", $energy,$kpoint,$band,$pdos[0],$ave,$diff_me,$pdos[6],$pdos[7],$diff_OCO;
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

printf "bands (min,max) = (%5d,%5d)\n", $band_min, $band_max;

sub average{
    my(@list)=@_;
    my($i,$sum);

    $sum=0;
    for($i=0;$i<=$#list;$i++){
	$sum+=$list[$i];
    }

    return ($sum/$i);
}

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


