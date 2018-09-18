#!/usr/bin/perl
# written by Joonho Park
# read a.xyz and write b.xyz
# translate the xyz coordinate w.r.t. an specified atom
# Usage: a.pl file.xyz 3[reference atom index starts from 1]

$fin=$ARGV[0];
$ref=$ARGV[1];

@fname=split(/\./,$fin);

if($fname[1] ne "xyz"){
    print "input error: the suffix should be xyz \n";
    exit(0);
}

$fout="$fname[0].orig.xyz";

open(IN,"<$fin");
open(OUT,">$fout");

@f_line=();
### read three atoms's coordinates
$i=0;
$j=0;

while($line=<IN>){
    push(@f_line,$line);
    if($i == $ref+1){
	@field=split(/\s+/,$line);
    	if($field[0] eq "") {shift(@field);}
	$x_shift=$field[1]; $y_shift=$field[2]; $z_shift=$field[3];
    }
}continue{
    $i++;
}
$max_line=$i;

for($i=0;$i<$max_line;$i++){
    if($i<2){ print $f_line[$i]; print OUT $f_line[$i]; next;}
    @field=split(/\s+/,$f_line[$i]);
    if($field[0] eq "") {shift(@field);}
    printf "$field[0]  %10.5f %10.5f %10.5f\n", $field[1]-$x_shift, $field[2]-$y_shift, $field[3]-$z_shift,"\n";
    printf OUT "$field[0]  %10.5f %10.5f %10.5f\n", $field[1]-$x_shift, $field[2]-$y_shift, $field[3]-$z_shift,"\n";
}

close(IN); 
close(OUT);
