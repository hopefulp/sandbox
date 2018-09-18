#!/usr/bin/perl
# written by Joonho Park
# read a.xyz and write b.xyz
# translate the xyz coordinate w.r.t. an specified atom
# Usage: a.pl file.xyz 1 12 13 [three reference atom indices starts from 1] debug
# debug for write process

use Math::Trig;

$fin=$ARGV[0];
$ref1=$ARGV[1];
$ref2=$ARGV[2];
$ref3=$ARGV[3];

for($i=0;$i<=$#ARGV;$i++){
    if($ARGV[$i] =~ /bug/){
	$DEBUG="YES";
    }
}

@fname=split(/\./,$fin);

@ref=();
if($#ARGV <0){
    print "test triangle is:\n";
#    @p1=qw(-1 0 0);
#    @p2=qw(0 1 0);
#    @p3=qw(0 0 1);
    @p1=qw(0.97 -8.81 3.19);
    @p2=qw(1.38 -7.4 4.58);
    @p3=qw(-0.74 -8.1 2.44);
}elsif($fname[1] ne "xyz"){
    print "input error: the suffix should be xyz \n";
    exit(0);
}

@ref=(\@p1,\@p2,\@p3);
#print join("  ",$\ref),"\n";
$fout="$fname[0].orig.xyz";

open(IN,"<$fin");
open(OUT,">$fout");

@f_line=();
### read three atoms's coordinates
$i=0;
$j=0;

### get three atom for plane
if($#ARGV >= 3){
while($line=<IN>){
    push(@f_line,$line);
    if($i<2){next;}
    @field=split(/\s+/,$line);
    if($field[0] eq "") {shift(@field);}
    if($i==$ref1+1 or $i==$ref2+1 or $i==$ref3+1){
	$ref[$j][0]=$field[1]; $ref[$j][1]=$field[2]; $ref[$j][2]=$field[3];
	$j++;
    }
    $kinds[$i-2]=$field[0]; $atoms[$i-2][0]=$field[1]; $atoms[$i-2][1]=$field[2]; $atoms[$i-2][2]=$field[3];
}continue{
    $i++;
}
$ntatom=$i-2;
}
if(defined $DEBUG){
    for($i=0;$i<3;$i++){
    	for($j=0;$j<3;$j++){
	    print $ref[$i][$j],"\t";
    	}
    	print "\n";
    }
}

for($i=0;$i<3;$i++){
#    print $ref[$i][0],$ref[$i][1],$ref[$i][2],"\n";
    $r12[$i]=$ref[1][$i]-$ref[0][$i];
    $r13[$i]=$ref[2][$i]-$ref[0][$i];
}
$nr[0]=$r12[1]*$r13[2]-$r12[2]*$r13[1];
$nr[1]=$r12[2]*$r13[0]-$r12[0]*$r13[2];
$nr[2]=$r12[0]*$r13[1]-$r12[1]*$r13[0];

$norm=0;
for($i=0;$i<3;$i++){
    $norm+=$nr[$i]*$nr[$i];
}
$norm=sqrt($norm);
#print $norm;
for($i=0;$i<3;$i++){
    $ei[$i]=$nr[$i]/$norm;
}

### ei test
#$eix=1/sqrt(3);
#@ei=(-$eix, $eix, $eix);
if(defined $DEBUG){    print "norm vector = ",join("  ",@ei),"\n";   }
### rotate phi w.r.t. y axis

#@pxy=qw($ei[0] $ei[1] 0);
$pxy[0]=$ei[0];
$pxy[1]=$ei[1];
$pxy[2]=0;

$yaxis[0]=0;
$yaxis[1]=1;
$yaxis[2]=0;
#print join("  ",@yaxis);
$xyinner=0;
for($i=0;$i<3;$i++){
    $xyinner+=$pxy[$i]*$yaxis[$i];    
}
$pxy_norm=sqrt($pxy[0]*$pxy[0]+$pxy[1]*$pxy[1]);
$yaxis_norm=abs($yaxis[1]);
$arg=$xyinner/($pxy_norm*$yaxis_norm);
#print "arg: $arg\t $xyinner\t $ei_norm\n";
#print "pi = ",pi,"   phi_arg = $phi_arg\n";
$phi=acos($arg);
#print "phi = ", $phi*360/pi, " \n";
if($pxy[0]<0){
    $phi=-$phi;
}
if(defined $DEBUG){ print "phi = $arg\t", $phi*360/pi, " \n"; }
#print cos($phi),"\t", acos(-1),"\n";

### rotation phi matrix
$cos11=cos($phi);
$sin12=-sin($phi);
$sin21=-$sin12;
@rot_phi1=($cos11, $sin12, 0);
@rot_phi2=($sin21, $cos11, 0);
@rot_phi3=qw(0 0 1);

#print join("  ",@rot_phi1),"\n";
@rot_phi=(\@rot_phi1,\@rot_phi2,\@rot_phi3);

for($i=0;$i<3;$i++){
    $ey[$i]=0;
    for($j=0;$j<3;$j++){
    	$ey[$i]+=$rot_phi[$i][$j]*$ei[$j];
    }
}
if(defined $DEBUG){ print "ey vector: ", join("  ",@ey),"\n"; }
### rotate all atoms for phi
for($n=0;$n<$ntatom;$n++){
    for($i=0;$i<3;$i++){
	$atatom_xy[$n][$i]=0;
	for($j=0;$j<3;$j++){
	    $atatom_xy[$n][$i]+=$rot_phi[$i][$j]*$atoms[$n][$j];
	}
    }
}
#for($n=0;$n<$ntatom+2;$n++){
#    if($n<2){ print $f_line[$n]; next;}
#    print $kinds[$n-2],"\t";
#    for($i=0;$i<3;$i++){
#	print  $atatom_xy[$n-2][$i],"\t";
#    }
#    print "\n";
#}
	


##### rotate theta w.r.t. x-axis to align with z-axis
$zaxis[0]=0;
$zaxis[1]=0;
$zaxis[2]=1;

$zinner=0;
for($i=0;$i<3;$i++){
    $zinner+=$ey[$i]*$zaxis[$i];
}
$zaxis_norm=abs($zaxis[2]);
$ey_norm=1;
$arg=$zinner/($zaxis_norm*$ey_norm);
$theta=acos($arg);
$cos22=cos($theta);
$sin23=-sin($theta);
$sin32=-$sin23;
if(defined $DEBUG) {print "theta: $theta\t",$theta*360/pi,"\n";}
@rot_theta1=qw(1 0 0);
@rot_theta2=(0, $cos22, $sin23);
@rot_theta3=(0, $sin32, $cos22);

#print join("  ",@rot_theta1),"\n";
@rot_theta=(\@rot_theta1,\@rot_theta2,\@rot_theta3);

for($i=0;$i<3;$i++){
    $ez[$i]=0;
    for($j=0;$j<3;$j++){
        $ez[$i]+=$rot_theta[$i][$j]*$ey[$j];
    }
}

if(defined $DEBUG) {print "ez vector: ", join("  ",@ez),"\n";}

### rotate all atoms for theta 
for($n=0;$n<$ntatom;$n++){
    for($i=0;$i<3;$i++){
	$atatom_new[$n][$i]=0;
	for($j=0;$j<3;$j++){
	    $atatom_new[$n][$i]+=$rot_theta[$i][$j]*$atatom_xy[$n][$j];
	}
    }
#    print $atatom_new[$n][0],$atatom_new[$n][1],$atatom_new[$n][2],"\n";
}
	
### write xyz
 
for($n=0;$n<$ntatom+2;$n++){
    if($n<2){ print $f_line[$n]; next;}
    print $kinds[$n-2],"\t";
    for($i=0;$i<3;$i++){
	print  $atatom_new[$n-2][$i],"\t";
    }
    print "\n";
}

close(IN); 
close(OUT);
