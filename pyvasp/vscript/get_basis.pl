#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

$atom=$ARGV[0];

if($#ARGV < 0){
    print "Error: input any atom symbol\n";
    exit(0);
}

$basfile="/home/joonho/basis/cc-pvtz.bas";
$outfile=$atom.".bas";
open(IN,"<$basfile");
open(OUT,">$outfile");
#print $basfile, $outfile;
$find="NO";
$flag="OFF";
@bas_line=();
### read three atoms's coordinates
#$i=0;
while($line=<IN>){
#    $i++;
    #print $line;
    @field=split(/\s+/,$line);
    if($field[0] eq "") {shift(@field);}
    if($flag eq "OFF"){			# wait until signal turns on
    	if($field[0] eq $atom){
#	    print $i."\n";
	    $flag="ON";
	    $find="YES";
	    push(@bas_line,$line);
	}else{ next;}
    }else{
 	if($line eq "****\n"){
	    push(@bas_line,$line);
	    $flag="OFF"; next;
	}else{
	    push(@bas_line,$line);}
    }
}

print @bas_line;
print OUT @bas_line;

if($find eq "NO"){
    print "Error: There is not atom $atom in $basfile\n";
}
close(IN); 
close(OUT);
