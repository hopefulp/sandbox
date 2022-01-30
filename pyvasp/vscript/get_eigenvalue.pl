#!/usr/bin/perl
# written by Joonho Park
# read Q-Chem output file, find key_word, print line
# Opt: after convergence
$fin=shift;
$nkey=shift || 1;
$nline=shift || 10;

open(IN,"<$fin");

$i=0;
$iline=0;
$ikey=0;
$j=0;
$After_word="OPTIMIZATION CONVERGED";
$tag_over="OFF";
$key_word="Alpha MOs";
$tag_count="OFF";
$ncount=0;
$vcount=0;
$stop_word="";

LINE: while($line=<IN>){
    $i++;
#    print $line;

    if($line !~ /$After_word/ and $tag_over eq "OFF"){
        next LINE;
    }elsif($line =~ /$After_word/){
	    print "over $After_word\n";
	    $tag_over="ON";         # this makes loop over this block after $tag_over="ON"
    }
	
    if($line =~ /$key_word/ and $tag_count eq "OFF"){
        $tag_count="ON";
        next LINE;
    }

    if($tag_count eq "ON"){
        $ncount++;
        if($ncount < 2){ print $line; next LINE; }
	    @line=split(/\s+/, $line);
        if($line[0] eq ""){ shift @line;}
        if( $line[1] =~ /^-\d/ or $line[1] =~ /^\d/){
            print "$ncount: $line";
        }elsif($line =~ /Virtual/){
            $vcount++;
            $tag_virtual="ON";
        }
        
        if($tag_virtual eq "ON" and $vcount < 3){
            $vcount++;
            next LINE;
        }elsif($tag_virtual eq "ON"){   # if protects block execution, last in case $tag_virtual is not "ON"
            last LINE;
        }
    }
}
continue{
    $i++;
}
close(IN);
    

