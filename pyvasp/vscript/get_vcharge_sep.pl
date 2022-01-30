#!/usr/bin/perl
# written by Joonho Park

use lib '/home/joonho/vaspi';
use MOF_constants;

$fin=shift(@ARGV);

@atom_list=();
while($#ARGV>=0){
    push(@atom_list,shift(@ARGV));
}

open(IN,"<$fin");

$i=0;
$iatom=1;
$sum=0;
while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}
    next if($line[0] =~ /\D/);
#    print $line[0]."\n";
    for($j=0;$j<=$#atom_list;$j++){
	if($line[0] == $atom_list[$j]){
	    $sum+=$line[4];
	    print "v_chg: $line[4]\n";
 	}
    }
} continue {
    $i++;
}
close(IN);

print $sum,"\n";


