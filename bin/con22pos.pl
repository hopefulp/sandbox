#!/usr/bin/perl
# written by Joonho Park
# CONTCAR (direct) POSCAR (cartesian) interconversion
# 

use Math::MatrixReal;

if($#ARGV < 0){
    print "Error: input CONTCAR or POSCAR";
    exit(0);
}

$fin=$ARGV[0];

@file=split(/\W+/,$fin);
$fname=$file[$#file];

if($fname eq "POSCAR"){
    $fout="CONTCAR";
}elsif($fname eq "CONTCAR"){
    $fout="POSCAR";
}else{
    print "error in input file\n";
    exit;
}

print "OUTFILE:: $fout\n";

open(IN, "<$fin");
open(OUT, ">$fout");
$i=0;
@atoms=();
@nsum=();
$ntatoms=0;
$iatom=0;
LINE: while($line=<IN>){
    ### copy 7 lines ###
    if($i < 7){
	print OUT $line;
	next LINE if($i < 1);
    }
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($i==1){ $scale=$line; next LINE;}
    if($i<5){
	for($j=0;$j<3;$j++){$lattice_vector[$i-2][$j]=$line[$j];}
	next LINE;
    }
    if($i==5){
	if($line[0] !~ /[A-Z]/) {
	    print "Error in $fin format: there are not atom symbols in 6th line\n";
	    exit(0);
	}else{ for($j=0;$j<=$#line;$j++){push(@atoms,$line[$j]);} }
	next LINE;
    }
    if($i==6){
	for($j=0;$j<=$#line;$j++){ 
	    $ntatoms+=$line[$j];
	    $nsum[$j]=$ntatoms;
 	}
	next LINE;
    }
    if($i==7){
	if($fname eq "POSCAR"){
	    die "It is not cartesian coordinate in $fname" unless($line[0] =~ /^[cC]/ );
	    print OUT "Direct\n";
	}elsif($fname eq "CONTCAR"){
    	    die "It is not direct coordinate in $fname" unless($line[0] =~ /^[dD]/ );
	    print OUT "Cartesian\n";
	}
    	next LINE;
    }
    
    for($j=0;$j<3;$j++){ $coord[$iatom][$j]=$line[$j];}
    $iatom++;
    if($iatom >= $ntatoms){ last;}		
} continue {
    $i++;
}
close(IN); 

$unit_matrix=Math::MatrixReal->new_from_cols([[$lattice_vector[0][0],$lattice_vector[0][1],$lattice_vector[0][2]],[$lattice_vector[1][0],$lattice_vector[1][1],$lattice_vector[1][2]],[$lattice_vector[2][0],$lattice_vector[2][1],$lattice_vector[2][2]]]);
$inv_matrix=$unit_matrix**-1;

if($file eq "CONTCAR"){
    for($i=0;$i<$ntatoms;$i++){
    	($x,$y,$z)=&xyz($coord[$i][0],$coord[$i][1],$coord[$i][2]);
#    	print $x*$scale,"\t",$y*$scale,"\t",$z*$scale,"\n";
    	print OUT $x*$scale,"\t",$y*$scale,"\t",$z*$scale,"\n";
#    	exit;
    }
}else{
    for($i=0;$i<$ntatoms;$i++){
	($ap,$bp,$cp)=&abc($coord[$i][0],$coord[$i][1],$coord[$i][2]);
	print OUT $ap,"\t",$bp,"\t",$cp,"\n";
    }
}

sub abc{
    my ($x, $y, $z)=@_;
    my ($ap,$bp,$cp,$x_col,$mat_mul);

    $x_col=Math::MatrixReal->new_from_cols([[$x,$y,$z]]);
    $mat_mul=$inv_matrix->multiply($x_col);

    $ap=$mat_mul->element(1,1);
    $bp=$mat_mul->element(2,1);
    $cp=$mat_mul->element(3,1);

    return ($ap,$bp,$cp);
}

sub xyz{
    my ($ap,$bp,$cp)=@_;
    my ($x, $y, $z);
    $x=$lattice_vector[0][0]*$ap+$lattice_vector[1][0]*$bp+$lattice_vector[2][0]*$cp;
    $y=$lattice_vector[0][1]*$ap+$lattice_vector[1][1]*$bp+$lattice_vector[2][1]*$cp;
    $z=$lattice_vector[0][2]*$ap+$lattice_vector[1][2]*$bp+$lattice_vector[2][2]*$cp;

    return ($x,$y,$z);
}
close(OUT);
