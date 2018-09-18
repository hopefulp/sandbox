#!/usr/bin/perl
# written by Joonho Park
# read  DOSCAR 

$fin=$ARGV[0];
#$sys=$ARGV[1];
$ref_sep=$ARGV[1];
#print $fin."\n";
$fout="ldos_peak.dat";

open(IN,"< $fin");
open(OUT,"> $fout");

#if($sys=~/(CO|co)2/){$ref_sep=1.e-3;}
#else {$ref_sep=0.0;}


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
for($i=0,$j=0;$i<$Nene;$i++){
#    printf "TDOS: Emin = %10.4f  Emax= %10.4f\n",$Emin, $Emax if($i==0);
    
    if($dos[$i] > $ref_sep){
	if($i!=0){
	    $dx=$energy[$i]-$energy[$i-1];
	}else{ $dx=$energy[$i+1]-$energy[$i];
	}
	$iene_idos=$energy[$i]*$dos[$i]*$dx;
	$idos=$dos[$i];
	$ave_num+=$iene_idos;
	$ave_den+=$idos*$dx;
	$j++;
    }elsif($ave_num!=0){
	print $ave_num/$ave_den,"\t",$ave_den,"\n";
	# resolve peak here!
	$ave_num=0; $ave_den=0;
    }

}

close(IN); close(OUT); #close(XYZ);
