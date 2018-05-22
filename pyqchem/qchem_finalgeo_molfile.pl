#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out
# write "mol.inp" as qchem molecule files

$fname=$ARGV[0];
@fname=split(/\./,$fname);
$fname_pre=$fname[0];
$fout="mol.inp";

$charge = 0;
$multip = 1;

if($#ARGV < 2){
    print "Warning: charge = $charge and multiplicity = $multip\n";
}

#print $#ARGV."\n";


######################### INPUT OF CHARGE AND SPIN (MULTIPICITY)

for($i=1;$i<=$#ARGV;$i++){
    @arg=split(/=/,$ARGV[$i]);
    if($arg[0] eq "charge"){
	$charge=$arg[1];
    }elsif($arg[0] eq "spin" or $arg[0] eq "multip"){
	$multip=$arg[1];
    }
}

#print $charge,$multip,"\n";
open(IN,"<$fname");
open(OUT,">$fout");

$nline=0;
$flag="OFF";
@line_list=();
$natom=100;
### read three atoms's coordinates
while($line=<IN>){
    @blank_test=split(/\s+/,$line);
    if($flag eq "ON" and $#blank_test == -1){ last;}
    if($flag eq "OFF"){
	#print $line;
    	@field=split(/\W+/,$line);
    	if($field[0] eq "") {shift(@field);}
    	if($field[0] ne "GEOMETRIES"){
	    next;
    	}else{ $flag="ON";}
    }elsif($flag eq "ON"){
	@field=split(/\s+/,$line);
	if($field[0] eq "") {shift(@field);}
	if($field[0] =~ /\d/){ 
	    if($nline%($natom+2) == 0){ @line_list=();}
	    $natom=$field[0]; 
	    #print $nline."\n";
	}
	push(@line_list,$line);	# save xyz file
	$nline++;
    }
}

print @line_list;
for($i=0;$i<=$#line_list;$i++){
    if($i<2) {next;}
    if($i==2) {print OUT "\$molecule\n$charge\t$multip\n";}
    print OUT $line_list[$i];
    if($i==$#line_list) {print OUT "\$end\n";}
}


