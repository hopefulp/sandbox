#!/usr/bin/perl
# written by Joonho Park
# read a.csh and write t.csh

$fpbs=$ARGV[0];
$fin=$ARGV[1];

#$group=$ARGV[2];
$ppn=$ARGV[2];

#if( $ppn eq ""){
#    if($group eq "g1"){ $ppn=12;}
#    elsif($group eq "g2"){ $ppn=16;}
#    else { print "group is not defined\n"; exit(1);}
#}

@base=split(/\./,$fin);
if($base[0] eq ""){ shift(@base);}
$basename=$base[0];

#@nodes=split(/=/,$nodes);
#$nnodes=$nodes[1];

$fout="t.csh";

open(IN,"<$fpbs");
open(OUT,">$fout");

@f_line=();
### read three atoms's coordinates
$i=0;
LINE: while($line=<IN>){
#    $i++;
    #print $line;
    if($i==0){
	print OUT "#PBS -N $basename\n";
	next LINE;
    }elsif($i==1){
#	if($nnodes eq ""){
#	    @part=split(/=/,$line);
#	    @nodes=split(/:/,$part[1]);
#	    $nnodes=$nodes[0];
#	}
#	if($group eq "g1"){
#	    
    	    print OUT "#PBS -l nodes=1:ppn=$ppn\n"; # :$group\n";
#	}else{
#	    print OUT "#PBS -l nodes=1:ppn=4:g2\n";
#	}
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
print "The job will be running on $nnodes nodes with $ppn processes \n"; # in $group group\n";

