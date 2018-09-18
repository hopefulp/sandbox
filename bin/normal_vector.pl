#!/usr/bin/perl
# written by Joonho Park
# read a.xyz and write b.xyz
# translate the xyz coordinate w.r.t. an specified atom
# Usage: a.pl file.xyz 3[reference atom index starts from 1]

$fin=$ARGV[0];
$ref1=$ARGV[1];
$ref2=$ARGV[2];
$ref3=$ARGV[3];

@fname=split(/\./,$fin);

@ref=();

if($#ARGV <0){
    print "test triangle is:\n";
    @p1=qw(-1 -1 0);
    @p2=qw(0 1 1);
    @p3=qw(1 0 -1);
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
if($#ARGV  => 3){
    while($line=<IN>){
        push(@f_line,$line);
        if($i==$ref1+1 or $i==$ref2+1 or $i==$ref3+1){
            @field=split(/\s+/,$line);
        	if($field[0] eq "") {shift(@field);}
            $ref[$j][0]=$field[1]; $ref[$j][1]=$field[2]; $ref[$j][2]=$field[3];
            $j++;
        }
    }continue{
    	$i++;
    }
    $max_line=$i;
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

### obtained plane vector (norm vector)
print "norm vector = ",join("  ",@ei),"\n";   

for($i=0;$i<$max_line;$i++){
    if($i<2){ print $f_line[$i]; print OUT $f_line[$i]; next;}
    @field=split(/\s+/,$f_line[$i]);
    if($field[0] eq "") {shift(@field);}
    printf "$field[0]  %10.5f %10.5f %10.5f\n", $field[1]-$x_shift, $field[2]-$y_shift, $field[3]-$z_shift,"\n";
    printf OUT "$field[0]  %10.5f %10.5f %10.5f\n", $field[1]-$x_shift, $field[2]-$y_shift, $field[3]-$z_shift,"\n";
}

close(IN); 
close(OUT);
