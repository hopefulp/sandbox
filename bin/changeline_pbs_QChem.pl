#!/usr/bin/perl
# written by Joonho Park
# read a.csh and write t.csh

$fin=$ARGV[0];
$infile=$ARGV[1];
$ppn=$ARGV[2];
#$group=$ARGV[3];

if($ppn eq ""){
    $ppn=12;
}

$nnodes=1;
$fout="t.csh";

open(IN,"<$fin");
open(OUT,">$fout");

@f_line=();
### read three atoms's coordinates
$i=0;
LINE: while($line=<IN>){
#    $i++;
    #print $line;
    if($i==0){
	print OUT "#PBS -N $infile\n";
	next LINE;
    }elsif($i==1){
#	if($group eq "g1"){
#    	    print OUT "#PBS -l nodes=$nnodes:ppn=12\n";
#    	    print OUT "#PBS -l nodes=$nnodes:ppn=12:g1\n";
#	}elsif($group eq "g2"){
#	    print OUT "#PBS -l nodes=$nnodes:ppn=16:g2\n";
#	}else{
	print OUT "#PBS -l nodes=$nnodes:ppn=$ppn\n";
	next LINE;
    }else{
	print OUT $line;
	next LINE;
    }
}
continue{
    $i++;
}
close(IN); 
close(OUT);

system("qsub $fout");
#if($group eq "g1"){
    print "The job will be running on $nnodes nodes with $ppn processes\n";
#}elsif($group eq "g2"){
#    print "The job will be running on $nnodes nodes with 16 processes\n";
#}else{
#    print "The job will be running on $nnodes nodes with 8 processes in idft\n";
#}

