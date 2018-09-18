#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

$fin=$ARGV[0];
@fname=split(/\./,$inp);
if($fname[1] ne "inp"){
    print "input error: the suffix should be inp of qchem input file\n";
    exit(0);
}

$fout="$fname[0].xyz";

open(IN,"<$fin");
open(OUT,">$fout");

$find="NO";
$flag="OFF";
$key="\$molecule";
$key2="\$end";
@f_line=();
### read three atoms's coordinates
$i=0;
while($line=<IN>){
#    $i++;
    #print $line;
    @field=split(/\s+/,$line);
    if($field[0] eq "") {shift(@field);}
    if($flag eq "OFF"){			# wait until signal turns on
    	if($field[0] eq $key){
#	    print $i."\n";
	    $flag="ON";
	}else{ next;}
    }else{
	$i++;
	if($field[0] eq $key2){
	    $flag="OFF"; next;
	}elsif($i<2){
	    next;
	}else{
	    push(@f_line,$line);}
    }
}

print OUT $#f_line+1, "\n\n", @f_line;

close(IN); 
close(OUT);
