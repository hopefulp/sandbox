#!/usr/bin/perl
# written by Joonho Park
# read POSCAR and write lattice constants, angles, and volume

use Math::Trig;
#use Math::T

if($#ARGV < 0){
    print "Error: input POSCAR format";
    exit(0);
}

$fpos=$ARGV[0];

open(IN,"<$fpos");

$key_lattice1="A3";
$key_lattice2="B3";
$key_lattice3="C3";
$key_atom="ACL";
$key_coord="XYZ";
$lattice_label=3;
$atom_label=3;
$xyz_label=3;
$natom=0;

@alattice=();
@blattice=();
@clattice=();


$find="NO";
$flag="OFF";
@bas_line=();
### read three atoms's coordinates

$nline=0;
while($line=<IN>){
    if($nline<2){
	$nline++;
	next;
    }

    @field=split(/\s+/,$line);
    if($field eq ""){
	$blank=shift(@field);
    }
    $i=$nline-2;
    for($j=0;$j<3;$j++){
	$coord[$i][$j]=$field[$j];
    }

    if($i==2) {last;}
    $nline++;
}

close(IN); 
@a=();
@b=();
@c=();

for($i=0;$i<3;$i++){
    push(@a,$coord[0][$i]);
    push(@b,$coord[1][$i]);
    push(@c,$coord[2][$i]);
}    

#print "@b\n";

#($x,$y,$z)=&vector_cross(@b,@c);
@bcrossc=&vector_cross(@b,@c);


#print "@bcrossc\n";
$vol = &vector_dot(@a,@bcrossc);
printf "volume                    = %10.2f\n", $vol;

for($i=0;$i<3;$i++){
    $asq+=$a[$i]*$a[$i];
    $bsq+=$b[$i]*$b[$i];
    $csq+=$c[$i]*$c[$i];
}

$arccos_a= vector_dot(@a,@b)/sqrt($asq*$bsq);
$arccos_b= vector_dot(@b,@c)/sqrt($bsq*$csq);
$arccos_c= vector_dot(@c,@a)/sqrt($csq*$asq);
#$cosa=atan2(sqrt(1.-$cosinva*$cosinva),$cosinva);
$alpha=acos($arccos_a);
$beta =acos($arccos_b);
$gamma=acos($arccos_c);
printf "lattice constants (a b c) = %10.4f   %10.4f   %10.4f\n", sqrt($asq), sqrt($bsq), sqrt($csq);
printf "angle  (alpha beta gamma) = %10.2f   %10.2f   %10.2f\n", $alpha*180/pi, $beta*180/pi, $gamma*180/pi ;




sub vector_dot{
    my (@a,@b);
    my ($i,$adotb);

    @a=();
    @b=();
    for($i=0;$i<3;$i++){
        push(@a,$_[$i]);
        push(@b,$_[$i+3]);
    }
    
    $adotb=$a[0]*$b[0]+$a[1]*$b[1]+$a[2]*$b[2];

    return $adotb;
}

sub vector_cross {
    my (@b,@c) ;
    my ($e1,$e2,$e3,$i);

    @b=();
    @c=();
    for($i=0;$i<3;$i++){
	push(@b,$_[$i]);
    	push(@c,$_[$i+3]);
    }

#    print "@b\n"; 
    $e1=$b[1]*$c[2]-$b[2]*$c[1];
    $e2=$b[2]*$c[0]-$b[0]*$c[2];
    $e3=$b[0]*$c[1]-$b[1]*$c[0];

    return ($e1,$e2,$e3);
}

