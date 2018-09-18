#!/usr/bin/perl
# written by Joonho Park
# read Q-Chem output file, find key_word, print line
# Opt: after convergence

use Data::Dumper qw(Dumper);
use Data::Dump qw(dump); 

$fin=shift;
### for 1-PP use (4, 7) and (4, 9) for CO2
### for 2-PPP (
$nLumo=4;       # for d10, Lumo is s, p    
$nHomo=10;       # for d10, Homo is d, ligand 2 for PP, 3 for PPP, PNP for sigma
#$ene_homo=shift;
$tag_print="NO";

open(IN,"<$fin");

$i=0;
$key_word="Alpha MOs";
$tag_count="OFF";
$ncount=0;
$vcount=0;

@eigen_values=();
$orbital="d";

LINE: while($line=<IN>){
    $i++;
#    print $line;

    if($line =~ /$key_word/ and $tag_count eq "OFF"){
        $tag_count="ON";
        next LINE;
    }

    if($tag_count eq "ON"){
        $ncount++;
        if($ncount < 2){ if($tag_print eq "YES") {print $line;}
            next LINE;
        }
	    @line=split(/\s+/, $line);
        if($line[0] eq ""){ shift @line;}
        #if( $line[0] =~ /^-\d/ or $line[1] =~ /^\d/){
        if($tag_print eq "YES") {print "$ncount: $line";}
        if($line !~ /Virtual/ and $vcount == 0){    # if $vcount != 0, skip this block
            for( $j=0 ; $j<=$#line; $j++){
                #if($line[$j] =~ /^\*/) { next; }
                if($tag_print eq "YES") {print $line[$j];}
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
                if($tag_print eq "YES") {print $line[$j];}
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

##### obtain LUMO energy
for($i=0; $i<=$#eigen_values; $i++){
    if($tag_print eq "YES") {printf "%d  %7.3f\n", $i+1, $eigen_values[$i];}
    if($eigen_values[$i] == $ene_lumo){ $l_lumo=$i; }
}

##### make hash
%LH_level=();
# $l_lumo level is lumo-1, +5(i=5) = lumo-6, 
for($i=$nLumo-1; $i>=-$nHomo; $i--){    # for i=-1 for Homo, 0 for Lumo 
    $j=$l_lumo+$i;                      # j start for l_lumo
    if($i >= 0){
        $lumo_num=$i+1;                 # whin i=0,  LUMO-1
        if($tag_print eq "YES") {printf "LUMO-%d: %7.3f\n", $lumo_num, $eigen_values[$j];}
        printf "%7.3f  ", $eigen_values[$j];
        $l_key="LUMO-"."$lumo_num";
        #print $l_key."\n";
        $LH_level{$l_key}=$eigen_values[$j];
    }else{
        $homo_num=abs($i);
        if($tag_print eq "YES") {printf "HOMO-%d: %7.3f\n", $homo_num, $eigen_values[$j];}
        printf "%7.3f  ", $eigen_values[$j];
        $l_key="HOMO-"."$homo_num";
        $LH_level{$l_key}=$eigen_values[$j];
    }
}
print "\n";
#print Dumper \%level;

#foreach my $keys (sort keys %LH_level){
#    #print $keys, "\t", $level{$keys}, "\n";
#    print $LH_level{$keys}, "  ";
#}
#print "\n";



