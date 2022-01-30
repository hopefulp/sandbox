#!/usr/bin/perl
# written by Joonho Park
# read Q-Chem output file, find key_word, print line
# Opt: after convergence

use Data::Dumper qw(Dumper);
use Data::Dump qw(dump); 
use List::Util 1.33 'any';

if($#ARGV < 0){
    print "Usage:: $0 qchem.out num_of_LUMO[4] num_of_HOMO[10] crit_coeff[0.2] check_atom=\"Ni P O\"\n";
    exit;
}    
$fin=shift;
### for 1-PP use (4, 7) and (4, 9) for CO2
### for 2-PPP (
$nLumo=shift || 4 ;       # for d10, Lumo is s, p    
$nHomo=shift // 10;       # for d10, Homo is d, ligand 2 for PP, 3 for PPP, PNP for sigma
$ntotal_orb=$nHomo+$nLumo;
$coeff_crit=shift || 0.2;
$check_atom="Ni P O C";
$check_atoms=shift || $check_atom;
@check_atoms=split(/\s+/, $check_atoms);
if($check_atoms eq "") {shift @check_atoms;}
print "check atoms:@check_atoms\n"; 

$print_level=0;
$put_only_energy=0;

print "num of lumo $nLumo, homo $nHomo \n";

open(IN,"<$fin");

$i=0;
$key_word="Alpha MOs";
$tag_keyword="OFF";
$tag_HL="ON";
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
@basis_coeff=();

$ib_coeff=0;
$tmp=0;
$nbasis=0;
$i_basis=0;
LINE: while($line=<IN>){
    #if($i==300){ exit(0); }
    ### up to key_word=Alpha MO, skip line
    ### after key_word, tag on and skip if-block
    ### if there are many key_word in outfile, add more tag here
    if($line =~ /$key_word/ and $tag_keyword eq "OFF"){
        if($print_level) {print $line;}
        $tag_keyword="ON";
        next LINE;
    } ### do not use else        
    ###                         up to where 
    if($tag_keyword eq "ON" and $tag_HL eq "ON"){
        if($line =~ /Occupied/){
            $i_occ++;
            if($i_occ != 0){ $vcount=0; }
            $tag_occ="ON";
            if($print_level) {print "here: $i_occ-th occ block\n";}
            next LINE;
        ### read occupied orbital energies
        }elsif($tag_occ eq "ON" and $tag_vir eq "OFF"){
	        @line=split(/\s+/, $line);
            if($line[0] eq ""){ shift @line;}
            ### 
            if($line !~ /Virtual/){
                for($j=0;$j<=$#line;$j++){push @{$eigen_values[$i_occ]}, $line[$j];}
                next LINE;
            ### if meet Virtual, vcount is raised
            ### above blick is skipped by vcount==0 and below block is skipped by tag_vir
            }elsif($line =~ /Virtual/){
                $tag_occ="OFF"; $tag_vir="ON";
                next LINE;
            }
        ### read virtual orbital energy
        }elsif($tag_occ eq "OFF" and $tag_vir="ON" ){
	        @line=split(/\s+/, $line);
            if($line[0] eq ""){ shift @line;}
            ### only one line is saved
            if($vcount==0){ 
                if($print_level) {print "here: $i_occ-th virtual block\n";}
                $ene_lumo[$i_occ]=$line[0]; 
                $index_lumo[$i_occ]=$#{$eigen_values[$i_occ]}+1; # index is for array
                #printf "lumo array index = %d: orbital id = %d\n", $index_lumo[$i_occ],$index_lumo[$i_occ]+1;
                printf "lumo orbital ene %7.3f, id = %d\n",$ene_lumo[$i_occ], $index_lumo[$i_occ]+1;
                for($j=0;$j<=$#line;$j++){push @{$eigen_values[$i_occ]}, $line[$j];}
                #print Dumper $eigen_values[$i_occ];
                $vcount++;
                next LINE;
            }
            if($line =~ /\d/){ next LINE; }
            elsif($line =~ /-------/) { $tag_HL = "OFF"; next LINE; }
            else{ $tag_vir="OFF"; next LINE; }
        }
    ##### Last execution block for wrap up
    ##### if protects block execution, last in case $tag_virtual is not "ON"
    }elsif($tag_HL eq "OFF" and $tag_keyword eq "ON"){
        ### print Homo-Lumo order
        for($j=0; $j<=$i_occ; $j++) {$index_min[$j]=$index_lumo[$j]-$nHomo;}
        for($j=0; $j<=$i_occ; $j++){ 
            #for($k=0; $k<$ntotal_orb; $k++){ print $eigen_values[$j][$index_min[$j]+$k],"\t";}
            for($k=$ntotal_orb-1;$k>=0;$k--){ print $eigen_values[$j][$index_min[$j]+$k],"\t";}
            print "\n";
        }
        $tag_keyword="OFF";     ### Alpha MO block done
        $n_occ=$i_occ+1;
        $i_occ=0;
        if($put_only_energy) { last LINE; }
        next LINE;
    }
    #########################################
    ##### BLOCK FOR GETTING COEFFICIENT just save #####
    #########################################
    if($line =~/$keyw_coeff/) { # this is true only once
        $tag_coeff="ON";
        next LINE;
    }elsif($tag_coeff eq "ON"){
        #chomp($line);
        #print "$i : $tag_cblock\n";
        ### check serial number line
        #if($line eq "\n") {last LINE; }       # after all the block there is blank line

        if($line !~ /\./){      # serial number check
            $ib_coeff++;        # start from 1
            @line=split(/\s+/, $line);
            if($line[0] eq ""){ shift @line;}
            #print $line[$#line],"\t",$index_mini[$i_occ]+1,"\n";
            if(($line[$#line] < $index_min[$i_occ]+1) or ($index_min[$i_occ]+$ntotal_orb+1 < $line[0])){
                $tag_bsave="NO";
            }else{ $tag_bsave="YES"; 
                print $line;
                push @{$basis_coeff[$i_occ][$ib_coeff-1]}, $line; 
            }
        }elsif($tag_bsave eq "NO"){
            next LINE;
        ##### search here: block save is YES and not \n, not serial number line
        }elsif($line ne "\n"){      
            push @{$basis_coeff[$i_occ][$ib_coeff-1]}, $line;   # ib_coeff scans all the blocks
            next LINE;
        }if($line eq "\n"){ 
            if($n_occ==1){
                last LINE;
            }elsif($n_occ==2){
                if($i_occ==0 )  { $i_occ++; $ib_coeff=0;}
                elsif($i_occ==1){ last LINE; }
            }
            
        }
    }        
}
continue{
    $i++;
}
close(IN);

@orbital_id=();
@col=();
@row=();
for($i=0;$i<$n_occ;$i++){
    $orbital_id[$i]=$index_min[$i]+1;
    $col=$orbital_id[$i]%6;
    $row=int ($orbital_id[$i]/6);   
    if($col!=0) {$col[$i]=$col; $row[$i]=($row++);}
    else        {$col[$i]=6;    $row[$i]=$row;}
    print  "$i-spin: orbital id $orbital_id[$i]:$row[$i] row $col[$i] column\n";
}    

print  "Alpha MO\n";
for($i=0;$i<$n_occ;$i++){                       # alpha beta index
    if($i==1) {print "Beta MO\n";}
    for($j=0;$j<=$#{$basis_coeff[$i]};$j++){     # coeff block index of 6 orbitals
        print  "this block: $i, $j\n"; 
        for($k=0;$k<=$#{$basis_coeff[$i][$j]};$k++){
            if($k<=1) {print $basis_coeff[$i][$j][$k]; next; }
            ### line decomposition
            @line=split(/\s+/,$basis_coeff[$i][$j][$k]);
            if($line[0] eq "") {shift @line;}
            if($line[2] ne "f"){ #print $line[2];
                for($l=3;$l<=$#line;$l++){
                    ### print with criterion of coeff_crit in case of check_atoms
                    if(abs($line[$l])>=$coeff_crit) { 
                        $is_atom="NO";
                        for($m=0; $m<=$#check_atoms; $m++){
                            if($line[1] eq $check_atoms[$m]) {$is_atom="YES"; last;}
                        }
                        if($is_atom eq "YES"){
                            ### for CO2 print all atoms
                            if($check_atoms[0] eq "C" or $check_atoms[1] eq "O"){ print  $basis_coeff[$i][$j][$k]; last;}
                            if($line[1] ne "C") { print $basis_coeff[$i][$j][$k]; last;}
                            elsif($line[2]==1)  { print $basis_coeff[$i][$j][$k]; last;}
                        }
                    }
                }
            }
        }   
    }
}    

