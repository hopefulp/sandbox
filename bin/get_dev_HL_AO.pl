#!/usr/bin/perl
# written by Joonho Park
# read Q-Chem output file, find key_word, print line
# Opt: after convergence

use Data::Dumper qw(Dumper);
use Data::Dump qw(dump); 

$fin=shift;
### for 1-PP use (4, 7) and (4, 9) for CO2
### for 2-PPP (
$nLumo=shift || 4 ;       # for d10, Lumo is s, p    
$nHomo=shift // 10;       # for d10, Homo is d, ligand 2 for PP, 3 for PPP, PNP for sigma
$ntotal_orb=$nHomo+$nLumo;
$coeff_crit=0.2;

$tag_print="NO";
$put_only_energy=0;

print "num of lumo $nLumo, homo $nHomo \n";

open(IN,"<$fin");

$i=0;
$key_word="Alpha MOs";
$tag_keyword="OFF";
$tag_HL="OFF";
$tag_occ="OFF";
$tag_vir="OFF";
$lcount=0;
$vcount=0;
$i_occ=-1;

@eigen_values=();   # 2D array
#@eigen_values_b=();
@index_lumo=();
@ene_lumo=();
@index_min=();

%get_HL=();
$orbital="d";
### for orbital coefficients
$keyw_coeff="MOLECULAR ORBITAL COEFFICIENTS";
$tag_coeff="OFF";
$tag_cblock="OFF";
$ib_coeff=0;
$tmp=0;
$nbasis=0;
$i_basis=0;
LINE: while($line=<IN>){
    ### before key_word=Alpha MO, skip line
    ### after key_word, tag on and skip if-block
    ### if there are many key_word in outfile, add more tag here
    if($line =~ /$key_word/ and $tag_keyword eq "OFF"){  
        print $line;
        $tag_keyword="ON";
        next LINE;
    }
    ### 
    if($tag_keyword eq "ON" and $tag_HL eq "OFF"){
        if($line =~ /Occupied/){
            $i_occ++;
            if($i_occ != 1){ $vcount=0; }
            $tag_occ="ON";
            print "here: $i_occ-th occ block\n";
            next LINE;
        ### read occupied orbital energies
        }elsif($tag_occ eq "ON" and $tag_vir eq "OFF"){
	        @line=split(/\s+/, $line);
            if($line[0] eq ""){ shift @line;}
            ### 
            if($line !~ /Virtual/){
                for( $j=0 ; $j<=$#line; $j++){
                    #if($line[$j] =~ /^\*/) { next; }
                    #if($tag_print eq "YES") {print $line[$j],"\t";}
                    push @{$eigen_values[$i_occ]}, $line[$j]; 
                    #else         { push (@eigen_values_b, $line[$j]); }
                }
                #if($tag_print eq "YES") { print "\n"; }
                next LINE;
            ### if meet Virtual, vcount is raised
            ### above blick is skipped by vcount==0 and below block is skipped by tag_vir
            }elsif($line =~ /Virtual/){
                $tag_vir="ON";
                $tag_occ="OFF";
                next LINE;
            }
        ### read virtual orbital energy
        }elsif($tag_occ eq "OFF" and $tag_vir eq "ON"){
            print "here: $i_occ-th virtual block\n";
	        @line=split(/\s+/, $line);
            if($line[0] eq ""){ shift @line;}
            if($vcount==0){ 
                $ene_lumo[$i_occ]=$line[0]; 
                $index_lumo[$i_occ]=$#{$eigen_values[$i_occ]}; # index is for array
                #print "index lumo = $index_lumo";
            }
            for( $j=0 ; $j<=$#line; $j++){
                #if($tag_print eq "tt") {print $line[$j];}
                push @{$eigen_values[$i_occ]}, $line[$j]; 
                #else          { push (@eigen_values_b, $line[$j]); }
            }
            print Dumper $eigen_values[0];
            $vcount++;
            if($line !~ /\w/ or $vcount > 0){
                $tag_vir = "OFF";
            }
            next LINE;
        ##### Last execution block for wrap up
        ##### if protects block execution, last in case $tag_virtual is not "ON"
        }else{
            #$index_lumo= &get_LH(\@eigen_values, $ene_lumo);     # index from 0
            for($j=0; $j<=$i_occ; $j++) {$index_min[$j]=$index_lumo[$j]-$nHomo; print $index_min[$j];}
            $tag_HL="ON";            # get out this block after getting HL
            #print Dumper \@eigen_values;
            for($j=0; $j<=$i_occ; $j++){ 
                for($k=0; $k<$ntotal_orb; $k++){ print $eigen_values[$j][$index_min[$j]+$k],"\t";}
                print "\n";
            }
            if($put_only_energy) { last LINE; }
            next LINE;
        }
    }
    #########################################
    ##### BLOCK FOR GETTING COEFFICIENT #####
    #########################################
    if($line =~/$keyw_coeff/) { # this is true only once
        $tag_coeff="ON";
        next LINE;
    }elsif($tag_coeff eq "ON"){
        #chomp($line);
        #print "$i : $tag_cblock\n";
        ### check serial number line
        if($line eq "\n") {last LINE; }       # after all the block there is blank line
        @line=split(/\s+/, $line);
        if($line[0] eq ""){ shift @line;}
        if($line !~ /\./){
            if($line[$#line] < $index_min){             $tag_bsave="NO";}
            elsif($index_min+$ntotal_orb < $line[0]) {  last LINE;}
            ### this is head line for save block
            else{ $tag_bsave="YES"; 
                print $line;
                #for($j=0; $j <= $#line; $j++){
                #    if($line[$j] == $index_min) { 
                #        $i_line=$j; 
                #        print $i_line; 
                #        last;
                #    }
                #}
            }
        }elsif($tag_bsave eq "NO"){
            next LINE;
        ##### search here
        }else{
            if($line =~ /eigenvalue/) {print $line; next LINE;}
            $pline=0;
            for($j=0; $j<6; $j++){
                if((($line[1] eq "Ni") or ($line[1] eq "P")) and ($line[$j+4] > $coeff_crit or $line[$j+4] < -$coeff_crit)){ 
                    $pline=1; last; 
                }
            }
            if($pline){ print $line; }
        }
    }        
}
continue{
    $i++;
}
close(IN);


##### make subroutine to get Homo Lumo list, hash and index of lumo-0
sub get_LH {
    my ($list_eigenvalue, $Ene_lumo, %hash_HL) = @_;

    ##### obtain LUMO energy
    my @list_eigenvalue=@{$list_eigenvalue};
    my %hash_HL=%{$hash_HL};
    for($i=0; $i<=$#list_eigenvalue; $i++){
        if($tag_print eq "YES") { printf "in sub: %d  %7.3f\n", $i+1, $eigen_values[$i];}
        if($eigen_values[$i] == $Ene_lumo){ $l_lumo=$i; }
    }
    return $l_lumo;
}


##### make hash
%LH_level=();
# $l_lumo level is lumo-1, +5(i=5) = lumo-6, 
#    for($i=$nLumo-1; $i>=-$nHomo; $i--){    # for i=-1 for Homo, 0 for Lumo 
#        $j=$l_lumo+$i;                      # j start for l_lumo
#        if($i >= 0){
#            $lumo_num=$i+1;                 # whin i=0,  LUMO-1
#            if($tag_print eq "YES") {printf "LUMO-%d: %7.3f\n", $lumo_num, $eigen_values[$j];}
#            printf "%7.3f  ", $eigen_values[$j];
##            $l_key="LUMO-"."$lumo_num";
 #           #print $l_key."\n";
 #           $LH_level{$l_key}=$eigen_values[$j];
 #       }else{
 #           $homo_num=abs($i);
 #           if($tag_print eq "YES") {printf "HOMO-%d: %7.3f\n", $homo_num, $eigen_values[$j];}
 #           printf "%7.3f  ", $eigen_values[$j];
 #           $l_key="HOMO-"."$homo_num";
  #          $LH_level{$l_key}=$eigen_values[$j];
  #      }
  #  }
  #  print "\n";
#print Dumper \%level;
  
#foreach my $keys (sort keys %LH_level){
#    #print $keys, "\t", $level{$keys}, "\n";
#    print $LH_level{$keys}, "  ";
#}
#print "\n";
