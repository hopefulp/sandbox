#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

$fname=$ARGV[0];
@fname=split(/\./,$fname);
$fname_pre=$fname[0];
$fout="$fname_pre.mol";
$fxyz="$fname_pre.xyz";
print "Overwrite mol file: $fout and $fxyz\n";
open(IN,"<$fname");
open(OUT,">$fout");
open(OUT1,">$fxyz");

if($ARGV[1] eq ""){
    $n_tag=1;
}else{
    $n_tag=$ARGV[1];
}
$tag_spin="molecule";
$flag_spin="OFF";
$charge_ini=10;
$charge=$charge_ini;

$tag_geometry="GEOMETRIES";	# for MOLDEN output
$i_tag=0;
$nline=0;
$flag="OFF";
@line_list=();
$natom=100;
$i=0;
### read three atoms's coordinates
while($line=<IN>){
    $i++;
    @blank_test=split(/\s+/,$line);
    ### extract charge and spin
    if($flag_spin eq "OFF" and $charge == $charge_ini){
	if($line =~ /\$$tag_spin/ ){ print "found \$$tag_spin\n";
	    $flag_spin = "ON";
	    next;
	}
    }elsif($flag_spin eq "ON"){
	@field=split(/\s+/,$line);
	if($field[0] eq "") {shift(@field);}
	$charge=$field[0]; $multiplicity=$field[1];
	$flag_spin = "OFF";
	print $charge, "***",$multiplicity,"\n";
    }
    ### extract optimized molecular geometry: end of geometries are blank line
    if($flag eq "ON" and $#blank_test == -1){ last;}
    if($flag eq "OFF"){
	#print $line;
    	@field=split(/\W+/,$line);
    	if($field[0] eq "") {shift(@field);}
    	if($field[0] eq "$tag_geometry"){
	    $i_tag++;
	    if($i_tag == $n_tag){ $flag="ON";}
	    else{   next;}
    	}else{	next;}
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
#=cut
}

print @line_list;
for($i=0;$i<=$#line_list;$i++){
    print OUT1 $line_list[$i];
}

for($i=0;$i<=$#line_list;$i++){
    if($i<2) {next;}
    if($i==2) {print OUT "\$molecule\n"; print OUT "  $charge   $multiplicity \n";}
    print OUT $line_list[$i];
    if($i==$#line_list) {print OUT "\$end\n";}
}
close(IN);
close(OUT);
close(OUT1);

