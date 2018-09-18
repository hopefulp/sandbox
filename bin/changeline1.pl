#!/usr/bin/perl
# written by Joonho Park
# read a.csh and write t.csh

$fin=$ARGV[0];
$dir=$ARGV[1];
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
	print OUT "#PBS -N $dir\n";
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

#system("qsub $fout");

