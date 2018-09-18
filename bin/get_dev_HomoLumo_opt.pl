#!/usr/bin/perl
# written by Joonho Park
# read Q-Chem output file, find key_word, print line
# Opt: after convergence

use Data::Dumper qw(Dumper);
use Data::Dump qw(dump); 

$fin=shift;
$nLumo=4;       # for d10, Lumo is s, p
$nHomo=10;       # for d10, Homo is d, ligand 2 for PP, 3 for PPP, PNP for sigma
#$ene_homo=shift;

open(IN,"<$fin");

$i=0;
$After_word="OPTIMIZATION CONVERGED";
$tag_over="OFF";
$key_word="Alpha MOs";
$tag_count="OFF";
$ncount=0;
$vcount=0;

@eigen_values=();
$orbital="d";

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
        #if( $line[0] =~ /^-\d/ or $line[1] =~ /^\d/){
        print "$ncount: $line";
        if($line !~ /Virtual/ and $vcount == 0){    # if $vcount != 0, skip this block
            for( $j=0 ; $j<=$#line; $j++){
                #if($line[$j] =~ /^\*/) { next; }
                print $line[$j];
                push (@eigen_values, $line[$j]);
            }
            next LINE;
        }elsif($line =~ /Virtual/){                 # when $vcount >= 1, this block is skipped due to $line =~/Virtual/
            $vcount++;
            $tag_virtual="ON";
            next LINE;
        }
        
        if($tag_virtual eq "ON" and $vcount < 3){
            if($vcount==1){ $ene_lumo=$line[0]; }
            for( $j=0 ; $j<=$#line; $j++){
                print $line[$j];
                push (@eigen_values, $line[$j]);
            }
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
    
#print Dumper \@eigen_values; 

for($i=0; $i<=$#eigen_values; $i++){
    printf "%d  %7.3f\n", $i+1, $eigen_values[$i];
    if($eigen_values[$i] == $ene_lumo){ $l_lumo=$i; }
}
%level=();
# $l_lumo level is lumo-1, +5(i=5) = lumo-6, 
for($i=$nLumo; $i>=-$nHomo; $i--){
    $j=$l_lumo+$i;
    if($j >= $l_lumo){
        $key_num=$i+1;
        printf "LUMO-%d: %7.3f\n", $key_num, $eigen_values[$j];
        $l_key="LUMO-"."$key_num";
        #print $l_key."\n";
        $level{$l_key}=$eigen_values[$j];
    }else{
        $key_num=abs($i);
        printf "HOMO-%d: %7.3f\n", $key_num, $eigen_values[$j];
        $l_key="HOMO-"."$key_num";
        $level{$l_key}=$eigen_values[$j];
    }
}

#print Dumper \%level;

foreach my $keys (sort keys %level){
    print $keys, "\t", $level{$keys}, "\n";
}







