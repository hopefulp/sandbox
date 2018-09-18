#!/usr/bin/perl
# written by Joonho Park
# energy overlap between $atom1 $atom1l orbital and $atom2 $atom2l orbital
# energy is given by center and delta E or Emin and Emax

use lib '/qcfs/joonho/modules';
use VaspDos;

if($#ARGV<0){
    print "Usage :: $0 [e=]-2.5[:E_max] \"atom list 3 9 etc\" \n";
    exit(0);
}

$switch_density="ldos";		# ldos for atom, pdos for angular momentum
$fin="PROCAR";

### analyze the second argument of energy for scanning energy region
$str_energy=$ARGV[0];
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

### analyze the third argument of atoms to be analyzed
@atoms=();
$str_atoms=$ARGV[1];
if($#ARGV == 1){
    @l_atoms=split(/\s+/,$str_atoms);
    if($l_atoms[0] eq ""){ shift @l_atoms;}
    for($i=0;$i<=$#l_atoms;$i++){
    	push(@atoms,$l_atoms[$i]);
    }
}else{
    print "running for all atoms\n";
    #print "Input atoms array in \"   \" ";
    #exit(2);
    #@atoms=qw(1 3);		# for CO2 DeG-B4C
    #@atoms=qw(5 6);		# for DeG-OBC
    #@atoms=qw(1 25); 	# for MOF and CO2
}
$natoms=$#atoms+1;
if( 0 <= $#atoms ){
    print "atom list: ($natoms atoms) ", join("  ",@atoms), "\n";
}else{
    print "atom list: all atoms \n";
}

##### 2. charge density for total density of atom (ldos), 0.02
#####                   for partial density of atom, less than 0.01
#@density=qw(0.00 0.00 0.00 0.00 0.00 0.00 0.00 ); 	# 0.02 0.02 ; 0.002 0.002 ;	0.003 0.003 
#@density=qw(0.01 0.01 ); 	# 0.02 0.02 ; 0.002 0.002 ;	0.003 0.003 
$density_ldos=0.02;    # check 1
$density_pdos=0.01;
$density_limit=$density_ldos;
if(! defined(@density)){
    @density=();
    for($i=0;$i<=$#atoms;$i++){
	push(@density, $density_limit);
    }
}

print " density criteria: ", join("  ",@density),"\n";

##### 1. if @atoms_l is commented, it calculates total density (ldos) of the atom
#@atoms_l=qw(2 1);
@atoms_l=();		# empty list is not to be defined
@atom_kind=();
if($switch_density eq "pdos"){
    for($i=0;$i<=$#atoms;$i++){
    	if($atoms[$i] <= 6){
            push(@atom_kind,"M");
            push(@atoms_l,2);        # for d-orbital
    	}else{
            push(@atom_kind,"O");
            push(@atoms_l,1);        # for p-orbital
    	}
    }
    print "angular momentum: ", join("  ",@atoms_l),"\n";
}
@pdos=();
@ldos=();
@dos=();
@m_orbital=(['s'],['py','pz','px'],['dxy','dyz','dz2','dxz','dx2']);


##### read PROCAR and analyze
open(IN,"<$fin");
$i=0;
$iatom=0;
$sum=0;
$band=0;
$ncheck_energy=0;
$skip_tag="OFF";	# do not skip now
$tot_coeff=0;
$tot_sq_coeff=0;
$n_kb=0;
print "read PROCAR\n";
print  "Energy\t\tKpoint\tband\t";
for($j=0;$j<=$#atoms;$j++){
    print" dos[$j]\t";
}
print "\n";
while($line=<IN>){
    #print "i=$i ... ;$line";
    next if($i<1);
    chomp($line);
    @line=split(/\W+/,$line);
    next if($#line < 0);

    ### get important parameters from 2nd line
    if($i==1){
	@line=split(/\W+/,$line);
	if($line[0] eq ""){shift(@line);}
	$tot_kpoints=$line[3];
	$weight=1./$tot_kpoints;
	$tot_bands=$line[6];
	$tot_atoms=$line[9];
	$max_skip=$tot_atoms+4;
#	print $tot_kpoints, $tot_bands, $tot_atoms, "\n";
	next;
   } 

    next if($skip_tag eq "ON");

    # read the following line
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}
	
    if($line[0] eq "k-point"){
	$kpoint=$line[1]; next;
    }elsif($line[0] eq "band"){
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
	next;
    }

#    if($line[0] eq "ion"){
#	print "calculation start\n"; next;
#    }

    ############### following #########################
    ### only for the selected energy
    ### for the specific atoms, check the density populated at the same energy

    # test for wavefunction coefficient
    #print $i,"\n";
#    if($#atoms < 0){
#	if($line[0] == 1) {print "\t\tline sum\tPROCAR   total   total_sq\n";}
#  	if(0<$line[0] and $line[0]<$tot_atoms){
#	    $sum_coeff=0; $sum_sq_coeff=0;
#	    for($j=0;$j<9;$j++){
#		$sum_coeff+=$line[$j+1];
#		$sum_sq_coeff+=$line[$j+1]*$line[$j+1];
#	    }
#
#	    $tot_coeff+=$sum_coeff;
#	    $tot_sq_coeff+=$sum_sq_coeff;
#	    #print "atom $line[0]:  $sum_coeff\t\t$line[10]  $tot_coeff\n";
#	    print "atom $line[0]:  $sum_sq_coeff\t\t$line[10]  $tot_coeff  $tot_sq_coeff\n";
#	}elsif($line[0] == $tot_atoms){
#	    exit(0);
#	}
#   }
#   next;

    for($j=0;$j<=$#atoms;$j++){
    	if($line[0] == $atoms[$j]){
##   	    if(defined(@atoms_l) and $ncheck_energy==1){	($nc,$nd)=&VaspDos::lm_quantum($atoms_l[$j]);}
	    # make column list for add as for each atom
	    @column_list=();
   	    if(defined(@atoms_l)){	@column_list=&VaspDos::lm_quantum($atoms_l[$j],$m,1); }
	    else{push(@column_list,9);} # this include only orbitals by "shift(@list); in sub sum()"
	    $nc=$column_list[0]; $nd=$column_list[$#column_list]; # column list should be consecutive
	    #print "$j\t$atoms_l[$j]\t$nc\t$nd\n";
	    if($switch_density eq "ldos"){
	        $ldos[$j]=$line[10];	# last (11-th) column is sum of all the other columns 
	 	$dos[$j]=$ldos[$j];
	    }else{
	    	$pdos[$j]=&sum($nc,$nd,@line);
		$dos[$j]=$pdos[$j];
	    }
	    #print join("  ",@column_list),"\n";
	    #print "index: ",$nc,"\t",$nd,"\n";
    	}
    }
    #### Have read one (k,band) set of atoms
    if($line[0] == $tot_atoms){
#    	print  $energy,"\t", $kpoint,"\t",$band,"\t",join("  ",@pdos),"\t", join("  ",@ldos),"\n";
	$tag_density="ON";
	for($j=0;$j<=$#atoms;$j++){
	    if($dos[$j] < $density[$j]){ $tag_density="OFF"; last;}
	}
	if($tag_density eq "ON"){
    	    print  "$energy\t$kpoint\t$band\t";
	    for($j=0;$j<=$#atoms;$j++){ print "$dos[$j]\t"; 	}
	    print "\n";
	    push @{ $kb_dos[$n_kb]}, @dos; 
	    $n_kb++;
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

### output summation
@kbsum_dos=();
for($i=0;$i<$n_kb;$i++){
#    print join("  ",@{ $kb_dos[$i]} ),"\n";
    for($j=0;$j<=$#atoms;$j++){
	$kbsum_dos[$j]+=$kb_dos[$i][$j];
    }
}

print "\t\t\t     ";
for($i=0;$i<=$#atoms;$i++){
    printf "%8.3f",$kbsum_dos[$i]*$weight;
}
print "\n";

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


