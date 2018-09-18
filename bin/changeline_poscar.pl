#!/usr/bin/perl
# written by Joonho Park
# read CONTCAR and write POSCAR

$fin=$ARGV[0];
$dir=$ARGV[1];
$fout="POSCAR";

open(IN,"<$fin");
open(OUT,">$fout");

@f_line=();
### read three atoms's coordinates
$i=0;
LINE: while($line=<IN>){
#    $i++;
    #print $line;
    if($line =~ /^[sS]/){
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

