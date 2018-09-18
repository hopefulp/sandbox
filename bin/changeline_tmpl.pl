#!/usr/bin/perl
# written by Joonho Park
# read incar and write t.incar
# usage:: ~ incar Me sys
# $ARGV[0] = incar.520.gga or incar.520.lda
# $ARGV[1] = Metal
# $ARGV[2] = cell or CO2; cell-MOF, CO2-MOF + 6CO2


$fin=$ARGV[0];
$Me=$ARGV[1];
$sys=$ARGV[2];

$fout="t.incar";

open(IN,"<$fin");
open(OUT,">$fout");

#@key_word=split(/\s+/,$line_arg);
$key_word="MAGMOM";

%mag_moment=(	"Mg" => 1.0, 	"Ca" => 1.0, 
		"Sc" => 1.5,	"Cu" => 1.5,
		"Ti" => 3,	"Ni" => 3,
		"V"  => 4.5,	"Co" => 4.5,
		"Cr" => 6,	"Fe" => 6,
		"Mn" => 7.5,
		"Zn" => 1.0
	    );

# for cell "48\*1" ; for CO2 "66\*1"
if($sys eq "cell"){
    $mag_later= "48\*1";
}elsif($sys eq "CO2"){
    $mag_later= "66\*1";
}elsif($sys eq "1CO2"){
    $mag_later= "51\*1";
}else{ print "there is not number of atoms signature\n"; exit();}

$mag_value=$mag_moment{$Me};
$new_line="$key_word = 6\*$mag_value $mag_later \n";

@f_line=();
### read three atoms's coordinates
$i=0;
LINE: while($line=<IN>){
#    $i++;
#    print $line;

    if($line =~ $key_word){
	$line=$new_line;
    }

    print OUT $line;    
    next LINE;
}
continue{
    $i++;
}
close(IN); 
close(OUT);

