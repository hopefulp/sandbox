#!/usr/bin/perl
# written by Joonho Park
# read a.xyz and write a0.xyz a1.xyz etc

$fname_pre=$ARGV[0];
$iatom[0]=$ARGV[1];	# number of line where an atom on the plane
$iatom[1]=$ARGV[2];	# number of line where an atom on the plane
$iatom[2]=$ARGV[3];	# number of line where an atom on the plane

$nAmol=$ARGV[4];
$nBmol=$ARGV[5];

$inscan=$ARGV[6];
$outscan=$ARGV[7];
$del=$ARGV[8];

for($i=0,$j=$inscan;$j<=$outscan;$i++,$j+=$del){
    }
$npoint=$i;
print "number of points = ",$npoint,"\n";

#$npoint=($outscan-$inscan)/$del;
#$npoint=$ARGV[8];
#$del=($outscan-$inscan)/($npoint-1);

$fname=$fname_pre.".xyz";
open(IN,"<$fname_pre");
open(OUT,">$fout");

$fA_a=$fname_pre."_a";
$fA_ab=$fname_pre."_ab";
$fB_b=$Bname."_b";
$fB_ab=$Bname."_ab";

$inp1=$fA_a.".inp";	$out1="c".$fA_a.".out";
$inp2=$fA_ab.".inp";	$out2="b".$fA_ab.".out";
$inp3=$fB_b.".inp";	$out3="e".$fB_b.".out";
$inp4=$fB_ab.".inp";	$out4="d".$fB_ab.".out";
$inp5=$fname;		$out5="a".$fname_pre.".out";

open(OUT0,">$inp1");
open(OUT1,">$inp2");
open(OUT2,">$inp3");
open(OUT3,">$inp4");

$i=0;
$nthree=0;
$tag=$i;
$flag="OFF";
$nskip=1;
@line_list=();
@x=();
@y=();
@z=();
### read three atoms's coordinates
while($line=<IN>){
    print "$i\n";
#    print "0:",$field[0],"1:",$field[1],"2:",$field[2],"\n";
    if($nthree >= 3) {last;}
    if($i==$iatom[0] || $i==$iatom[1] || $i==$iatom[2]){
	print $line;
    	@field=split(/\s+/,$line);
	if($field[0] eq "") {shift(@field);}
	push(@x,$field[1]);
	push(@y,$field[2]);
	push(@z,$field[3]);
	$nthree++;
        $i++;
    }else{$i++; next;}
    
}

### get plane vector
$A=$y[0]*($z[1]-$z[2])+$y[1]*($z[2]-$z[0])+$y[2]*($z[0]-$z[1]);
$B=$z[0]*($x[1]-$x[2])+$z[1]*($x[2]-$x[0])+$z[2]*($x[0]-$x[1]);
$C=$x[0]*($y[1]-$y[2])+$x[1]*($y[2]-$y[0])+$x[2]*($y[0]-$y[1]);
$mag=sqrt($A*$A+$B*$B+$C*$C);
print "A=",$A,"\tB=", $B,"\tC=", $C,"\tMagnitude=", $mag,"\n";
if($mag <= 0.0){
    print "error: directional vector has zero magnitude\n";
    exit(0);
}
$A/=$mag;
$B/=$mag;
$C/=$mag;

print "directional vector of plane", $A, $B, $C, "\n";

###
close(IN);
#close(OUT0);	close(OUT1);	close(OUT2); 	close(OUT3);
#system("qchem $inp1 $out1 &");
#system("qchem $inp2 $out2 &");
#system("qchem $inp3 $out3 &");
#system("qchem $inp4 $out4 &");
#system("qchem $inp5 $out5 &");

sub print_file {
    print OUT0 $line;
    print OUT1 $line;
    print OUT2 $line;
    print OUT3 $line;
}
