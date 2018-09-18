#!/usr/bin/perl
# written by Joonho Park
# read a.csh and write t.csh

$fin=$ARGV[0];
$infile=$ARGV[1];
$ppn=$ARGV[2];

if($ppn eq ""){
    $ppn=16;
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
print "The job will be running on $nnodes nodes with $ppn processes\n";

