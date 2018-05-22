#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

$fname=$ARGV[0];
@fname=split(/\./,$fname);
$fname_pre=$fname[0];
#$fout="mol.inp";

open(IN,"<$fname");
#open(OUT,">$fout");

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
#for($i=0;$i<=$#line_list;$i++){
#    if($i<2) {next;}
#    if($i==2) {print OUT "\$molecule\n";}
#    print OUT $line_list[$i];
#    if($i==$#line_list) {print OUT "\$end\n";}
#}


