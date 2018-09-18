#!/usr/bin/perl
# written by Joonho Park
# read LDOS file 
# usage ~ dos.file 1.e-5 Me [ads|des]
# ref-sep to seperate the peak
# Me and [ads|des] for different Fermi level.

use lib '/home/joonho/vaspi';
use MOF_constants;


$fin=$ARGV[0];
#$sys=$ARGV[1];
$ref_sep=$ARGV[1];
#print $fin."\n";
$fout="ldos_peak.dat";

$Me=$ARGV[2];
$sorp=$ARGV[3];

$fermi=&MOF_constants::get_fermi("MOF-1co2",$Me,$sorp);
#print $fermi,"\n";
#exit;

open(IN,"< $fin");
open(OUT,"> $fout");

$i=0;
@energy=();
@dos=();
### read three atoms's coordinates
LINE: while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}
    push(@energy,$line[0]);
    push(@dos,$line[1]);

### DOS LINE: 1. energy 2. tdos_@_E    3. tdos_int_E 						for non-spin
### DOS LINE: 1. energy 2. tdos_@_E_up 3. tdos_@_E_down 4. tdos_int_E_up 5. tdos_int_E_down 	for spin

} continue {
    $i++;
}
$Nene=$i;
# other options
# $ene_Fermi = 0;	# do not shift w.r.t Fermi level

$p_ene=0;
$i_peak=0;
$i_dx=0;
$t_dos_num=0.0;
$t_dos_den=0.0;
$t_anti_num=0.0;
$t_anti_den=0.0;
$tot_ene=0.0;
$dx=($energy[$#energy]-$energy[0])/($#energy-0);
#print "dx= $dx\n";
for($i=0,$j=0;$i<$Nene;$i++){
#    printf "TDOS: Emin = %10.4f  Emax= %10.4f\n",$Emin, $Emax if($i==0);

    $iene_idos=$energy[$i]*$dos[$i];
    $idos=$dos[$i];
    if($energy[$i]<=$fermi){
	$t_dos_num+=$iene_idos;
	$t_dos_den+=$idos;	
    }else{
	$t_anti_num+=$iene_idos;
	$t_anti_den+=$idos;	
    }
    ### for peak analysis    
    if($dos[$i] > $ref_sep){
#	$iene_idos=$energy[$i]*$dos[$i]*$dx;
#	$idos=$dos[$i]*$dx;
	$ave_num+=$iene_idos;
	$ave_den+=$idos;
	$j++;
    }elsif($ave_num!=0){
#	print $ave_num/$ave_den,"\t",$ave_den*$dx,"\n";
	print $i_peak+1,"\t",$ave_num/$ave_den,"\t",$ave_den*$dx,"\n";
	# resolve peak here!
	$ave_num=0; $ave_den=0; $i_peak++;
    }

}
print "Total energy = ",$t_dos_num*$dx,"\n";
print "Sum to Fermi: ",$t_dos_num/$t_dos_den,"\t",$t_dos_den*$dx,"\t";
print "Sum upper Fermi: ",$t_anti_num/$t_anti_den,"\t",$t_anti_den*$dx,"\n";
close(IN); close(OUT); #close(XYZ);
