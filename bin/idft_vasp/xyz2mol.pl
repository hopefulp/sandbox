#!/usr/bin/perl
# written by Joonho Park
# read a.mol and write a.xyz

$fin=shift @ARGV;
$suffix="mol";

$charge=0;
$multi=1;

use lib '/home/joonho/vaspi';
use IOFile;
$fout=file1($fin,$suffix);

if(1<=$#ARGV){ $charge=$ARGV[1];}
if(2<=$#ARGV){ $multi=$ARGV[2];}

open(IN,"<$fin");
open(OUT,">$fout");

$find="NO";
$flag="OFF";

print OUT "\$molecule\n";
print OUT "\t$charge\t$multi\n";

@f_line=();
### read three atoms's coordinates
$i=0;
while($line=<IN>){
    $i++;
    if($i <= 2){next;} 
    else{    push(@f_line,$line);}
    
}

print OUT @f_line;
print OUT "\$end";

close(IN); 
close(OUT);
