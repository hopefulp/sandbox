#!/usr/bin/perl
# written by Joonho Park
# read a.csh and write t.csh

if($#ARGV < 1){
    print "Usage:: $0 pbs-filename jobname\n";
    exit(1);
}

$fin=$ARGV[0];
$jobname=$ARGV[1];

$fout="t.csh";

open(IN,"<$fin");
open(OUT,">$fout");

$i=0;
LINE: while($line=<IN>){
    @line=split(/\W+/,$line);

    #if($i < 5) {print $line[2],"\n";}
    if($line[2] eq "N"){
        #print "$i-line:  $line";
	print OUT "#PBS -N $jobname\n";
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

