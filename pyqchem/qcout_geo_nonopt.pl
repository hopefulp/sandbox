#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

$fname=$ARGV[0];
@fname=split(/\./,$fname);
$fname_pre=$fname[0];
$fout="$fname_pre.non.mol";

open(IN,"<$fname");
open(OUT,">$fout");

$nline=0;
$find_flag="OFF";
@line_list=();
$natom=100;
$nmargin=3;
$iopt=0;
### read three atoms's coordinates
FLINE:while($line=<IN>){
    @geo_line=split(/\s+/,$line);
    @field=split(/\W+/,$line);
    if($field[0] eq "") {shift(@field);}
    #if($find_flag eq "ON" and $#geo_line == -1){ last;}
    if($find_flag eq "OFF"){
	#print $line;
    	if($field[0] ne "Optimization"){
	    next FLINE;			#skip before block
    	}else{ 
	    $find_flag="ON";		#start block
	    @line_list=();		#empty atom list
	    $iopt_o=$iopt;
	    $iopt=$field[2];
	    if($iopt < $iopt_o){ last;} #if another Opt process, exit
	}
    }elsif($find_flag eq "ON"){
	if($field[0] !~ /\d/){		#skip in the block
	    if($start_flag eq "ON"){		#exit block when copy done
	 	$find_flag="OFF";
		$start_flag="OFF";	#initialize start flag in the block
	    }
	    next FLINE;
	}else{
	    $start_flag="ON";
	    push(@line_list,$line);
	}
	$nline++;
    }
}
printf "Otimization step %d\n", $iopt-1;
#print @line_list;
WRITE:for($i=0;$i<=$#line_list;$i++){
    @field=split(/\s+/,$line_list[$i]);
    if($field[0] eq "") {shift(@field);} 
    if($i==0) {
	printf "%d\n\n", $#line_list+1; 
	print OUT "\$molecule\n"; print OUT "  0   1 \n";
    }	
    printf OUT "%2s\t%10.5f\t%10.5f\t%10.5f\n", $field[1], $field[2], $field[3], $field[4];
    printf "%2s\t%10.5f\t%10.5f\t%10.5f\n", $field[1], $field[2], $field[3], $field[4];
    if($i==$#line_list) {print OUT "\$end\n";}
}


