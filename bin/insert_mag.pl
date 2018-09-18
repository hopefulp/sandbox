#!/usr/bin/perl
# written by Joonho Park
# read incar and write t.incar
# usage:: ~ incar Me sys
# $ARGV[0] = INCAR w or w/o MAGMOM
# $ARGV[1] = Metal
# $ARGV[2] = number of Metal
# $ARGV[3] = system or # of all the other atoms
# $ARGV[4] = magnetic order
if($#ARGV < 3){
    print "Usage:: $0 INCAR_file Metal nmetal Magnetism order\n";
    print "   Magnetism:: FM for overall FM\n";
    print "               AFM for FM/AFM; inner-chain FM and inter-chain AFM\n";
    print "               IAFM for inner-chain AFM\n";
    print "               read for reading dat file\n";
    exit(1);
}

$fin=$ARGV[0];
$Me=$ARGV[1];
$nmetal=$ARGV[2];
$Mag=$ARGV[3];
$mag_order=$ARGV[4];	# double

$Me_array="half";	# alternative, half
$Me_SC="";	# double 1 2 3 7 8 9 	4 5 6 10 11 12

$Mag="AFM" if ! defined $Mag; 
if($Mag =~ /^i/i){
    if($mag_order eq ""){
	print "There is no input for magnetic ordering for Internal AFM\n";
	exit(1);
    }else{
	@mag_order=split(/\s+/,$mag_order);
	if($#mag_order >= 1){
	    print "Read magnetic ordering in command line\n";
	}else{
	    open(IN,"<$mag_order") or die "can't open $mag_order : $!";
	    $all_line="";
	    while($line=<IN>){
		$all_line.=$line;
	    }
	    close(IN);
	    @mag_order=split(/\s+/,$all_line);
	}
    }
}

print join("  ",@mag_order)."\n";

$fout="t.incar";

open(IN,"<$fin");
open(OUT,">$fout");

$key_magmom="MAGMOM";
$key_icharg="ICHARG";
$key_ispin="ISPIN";
$key_istart="ISTART";
%mag_moment=(	"Mg" => 1.0, 	"Ca" => 1.0, 
		"Sc" => 1.5,	"Cu" => 1.5,
		"Ti" => 3,	"Ni" => 3,
		"V"  => 4.5,	"Co" => 4.5,
		"Cr" => 6,	"Fe" => 6,
		"Mn" => 7.5,
		"Zn" => 1.0
	    );

$mag_value=$mag_moment{$Me};
$mag_odds="999*1.5";
$mag_odds0="999*0";
$mag_list="";
$nmetal_2=$nmetal/2;
if($Mag =~ /^F/i){
    $mag_list="$nmetal\*$mag_value $mag_odds";
}elsif($Mag =~ /^a/i){
    $mag_list="$nmetal_2\*$mag_value $nmetal_2\*-$mag_value ";
    if ($mag_order eq "double"){
	$mag_list.="$nmetal_2\*$mag_value $nmetal_2\*-$mag_value ";
    }
    $mag_list.="$mag_odds0 ";
}elsif($Mag =~ /^i/i){
    #### Alternative of up and down
    if($#mag_order != $nmetal-1){
	print "number of atoms and magnetic moments are not matched\n";
	exit (3);
    }
    for($i=0;$i<$nmetal;$i++){
	if($mag_order[$i] > 0){
	    	$mag_list.=" $mag_value ";
	}else{	$mag_list.=" -$mag_value ";
        }
    }  
    $mag_list.="$mag_odds0 ";
}else{
    print "Error: Magnetism is not defined\n";
    exit(10);
}

$new_line="$key_magmom = $mag_list \n";
print $new_line;

$i=0;
$tag_mod="OFF";
LINE: while($line=<IN>){

    if($line =~ $key_magmom){
	$line=$new_line;
	$tag_mod="ON";
    }elsif($line =~ $key_istart){
	$line="ISTART = 0\n";
    }elsif($line =~ $key_icharg){
	$line="ICHARG = 2\n";
    }elsif($line =~ $key_ispin){
	$line="ISPIN = 2\n";
    }

    if($line =~ "cont" and $tag_mod eq "OFF"){
	print OUT $new_line,"\n";
	$tag_mod="ON";
    }
	
    next LINE;
}
continue{
    $i++;
    print OUT $line;    
}
close(IN); 
close(OUT);
if($tag_mod eq "OFF"){
    print "Error:: MAGMOM was not modified\n";
    exit(10);
}

