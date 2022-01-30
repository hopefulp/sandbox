#!/usr/bin/perl
# written by Joonho Park

use lib '/home/joonho/vaspi';
use MOF_constants;

$fin=shift(@ARGV);

@arg_list=();
while($#ARGV>=0){
    push(@arg_list,shift(@ARGV));
}
@atom_index=();
@atom_list=();
while($#arg_list>=0){
    if($arg_list[0] =~ /\d/){
	$index=shift(@arg_list);
	push(@atom_index,$index);
    }else{
	$atom=shift(@arg_list);
	push(@atom_list,$atom);
    }
}
#print join(" ",@atom_list),"\n";

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
    for($j=0;$j<=$#atom_index;$j++){
	if($line[0] == $atom_index[$j]){
	    $v_chg=$line[4];
	    $sum+=$v_chg;
	    $p_chg=&MOF_constants::get_pristine_chg($atom_list[$j]);
	    printf "%10.4f \t", $p_chg-$v_chg;
#	    print $p_chg,"\n";
#	    print "v_chg: $line[4]\n";
 	}
    }
} continue {
    $i++;
}
print "\n";
close(IN);

#print $sum,"\n";


