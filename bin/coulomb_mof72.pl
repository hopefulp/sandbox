#!/usr/bin/perl
# CONTCAR.xyz is xyz file with the order of POSCAR

$fin=$ARGV[0];
open(IN,"<$fin");

$i=0;
$iatom=0;
$atoms=();
#### read all coordinates
while($line=<IN>){
    chomp($line);
    @line=split(/\s+/, $line);
    if($line[0] eq ""){ shift(@line);}
    if($i==0 and $line[0] =~ /\D/){
	print "Error: In xyz format, first line is expected to be a number!\n";
	exit;
    }
    next if($i<2);
     
    ### coords index starts from 0
    if($#line == 3){
        for($j=0;$j<3;$j++){
            $coords[$iatom][$j]=$line[$j+1];
        }
	$atoms[$iatom]=$line[0];
    	$iatom++;
    }else{
 	next;
    }
} continue {
    $i++;
}
$ntatom=$iatom;
close(IN);

#for($i=0;$i<$ntatom;$i++){
#    for($j=0;$j<3;$j++){
#    	print $coords[$i][$j],"\t";
#    }
#    print "\n";
#}

### second file is for Charges
$fin=$ARGV[1];
open(IN,"<$fin");

$i=0;
$iatom=0;
@charges=(); 	### charge indext start from 0
while($line=<IN>){
    chomp($line);
    @line=split(/\s+/, $line);
    if($line[0] eq ""){ shift(@line);}
    next if($i<2);

    if($line[0]-1 == $iatom){
	push(@charges,$line[4]);
	$iatom++;
    }else{ 
	next;
    }
}continue{
    $i++;
} 
close(IN);

if($iatom != $ntatom){
    print "Error in the number of atoms in $ARGV[0] $ntatom and $ARGV[1] $iatom\n";
    exit;
}

%pair_O1=( 1 => '32', 2 => '29', 3 => '33', 4 => '25', 5 => '27', 6 => '36');
%pair_O2=( 1 => '31', 2 => '30', 3 => '34', 4 => '26', 5 => '28', 6 => '35');
%pair_C=( 1 => 64, 2 => 63, 3 => 65, 4 => 61, 5 => 62, 6 => 66);
%ref_charge=( "O" => 6.0, "C" => 4.0, "H" => 1.0,"Mg"=>2.0,"Ca" => 10.0,"Sc"=>11.0,"Ti"=>4.0,"V"=>5.0,
		"Cr" => 6.0,"Mn"=>7.0,"Fe"=>8.0,"Co" => 9.0,"Ni"=>10.0,"Cu"=>11.0,"Zn"=>12.0  );
use lib '/home/joonho/modules';
use Arith;

$pair_size=0;
$pair_size+=scalar keys %pair_O1;
#print $pair_size."\n";

$coulomb_coeff=8.987551787368E9;		# N m^2 / C^2
$e_chg= 1.6021766E-19;				# C
$r_unit=1E-10;					# m
$coeff1=$coulomb_coeff*$e_chg*$e_chg/$r_unit;	# N m^2 /C^2 C^2 m^-1 = J (A)
$coeff=$coeff1*6.022141E23*1.E-3;
#printf "  %g  \n", $coeff;

$dist_sum=0;
$ECoulomb_sum=0;
### Metal index
for($i=0;$i<$pair_size;$i++){
    $tchg_me=$charges[$i];
    $chg_me=$ref_charge{$ARGV[2]}-$tchg_me;
	$ECoulomb=0;
    	#### pair index
	# 1st pair
	$p1=$i+1;	# particle index
	$j=$pair_O1{$p1}; #1print "Me = $p1, O1 = $j\n";
	$dist=&Arith::dist($coords[$i][0],$coords[$i][1],$coords[$i][2],$coords[$j-1][0],$coords[$j-1][1],$coords[$j-1][2]);
	error_dist($dist, $ARGV[2], $i+1,"O", $j);
	$ECoulomb+=get_energy($charges[$j-1], $ref_charge{"O"}, $coeff, $chg_me, $dist);
#  	print $dist,"\n";
	$dist_sum+=$dist;

	$j=$pair_O2{$p1}; #1print "Me = $p1, O1 = $j\n";
	$dist=&Arith::dist($coords[$i][0],$coords[$i][1],$coords[$i][2],$coords[$j-1][0],$coords[$j-1][1],$coords[$j-1][2]);
	error_dist($dist, $ARGV[2], $i+1, "O", $j);
	$ECoulomb+=get_energy($charges[$j-1], $ref_charge{"O"}, $coeff, $chg_me, $dist);

	$j=$pair_C{$p1};
	$dist=&Arith::dist($coords[$i][0],$coords[$i][1],$coords[$i][2],$coords[$j-1][0],$coords[$j-1][1],$coords[$j-1][2]);
	error_dist($dist, $ARGV[2], $i+1,"C", $j);
	$ECoulomb+=get_energy($charges[$j-1], $ref_charge{"C"}, $coeff, $chg_me, $dist);
#1	print "$chg_me $net_co2 $dist $ECoulomb\n";

#1	print "Coulomb energy between Me and CO2 $ECoulomb\n";
	$ECoulomb_sum+=$ECoulomb;

#	print "$coeff $chg_me $net_co2 $dist $ECoulomb\n";
#	$dist_sum+=$dist;
#	print " $tchg_me $chg_me $tchg_co2 $net_co2 \n";
#	printf "pair %d and %d\n",$i+1,$j;
}

#print "Ave distance is",  $dist_sum/6, "\n";
print $ECoulomb_sum/6, "\n";

sub get_energy{
    my($tchg_co2,$ref_charge,$pre_coeff,$chg1, $dist)=@_;
    my($net_co2,$net_energy);
    
    $net_co2=$ref_charge-$tchg_co2;
    $net_energy=$pre_coeff*$chg1*$net_co2/$dist;

#1    print "$chg1 $net_co2 $dist $net_energy\n";
#    print "$tchg_co2  $ref_charge  $pre_coeff $chg1 $chg2 $net_energy\n";
    return $net_energy;
}

sub error_dist{
    my($dist,$Me,$i,$atom,$j)=@_;
    my($crit);

    $crit=5.5;

    if($dist >= $crit){
	print "Error in measure distance between $Me $i and $atom $j as distance $dist\n";
	exit(0);
    }
}
