#!/usr/bin/perl
## written by Joonho Park
## read a.Amol, write a.mol a_ext_charges.mol

if($#ARGV < 0){
    print "Usage: $0 Amol-file charge charged_atom \n";
    exit (1);
}

$fin=shift @ARGV;
$ext_charge=shift @ARGV;
$ch_atom=shift @ARGV;
#print "atom list: ", join("  ",@list_atom) ,"\n";

$charge=0;
$multi=1;

use lib '/qcfs/joonho/modules';
use IOFile2;

# output file
($fout,$suffix)=IOFile2::file3($fin);
$fout1=$fout.".mol";
$fout2=$fout."_ext_charges.mol";

#if(0<=$#ARGV){ $charge=shift @ARGV;}
#if(0<=$#ARGV){ $multi=shift @ARGV;}


open(IN,"<$fin");
open(OUT1,">$fout1");
open(OUT2,">$fout2");

#print "fin: $fin; fout: $fout\n";

### read Accelrys mol file
@f_line=();
$i=0;
$natom=0;
$n_chatom=0;
@ch_coords=();
LOOP1: while($line=<IN>){
    # skip three lines
    if($i <= 2){next LOOP1;}
    else{
	@ncol=split(/\s+/,$line);
	if($ncol[0] eq ""){ shift @ncol;}
	if($i==3){ $nt_atom=$ncol[0]; next LOOP1;}
	if($i==4){ # print "\n\$molecule\n\t0\t1\n"; 
	    print OUT1  "\n\$molecule\n\t0\t1\n";
	    print OUT2  "\n\$molecule\n\t0\t1\n";
	}
#    	$atom_tag="N";
	$ch_atom_tag="N";
	#LOOP2: for($j=0;$j<=$#list_atom;$j++){
	if($ncol[3] =~ /\d/) { print "end of atom list\n"; last LOOP1;}
	if($ncol[3] eq $ch_atom){
		$ch_atom_tag="Y";
	    }
	}
	if($ch_atom_tag eq "N"){
	    #print  "$ncol[3]\t$ncol[0]\t$ncol[1]\t$ncol[2]\n";
	    print OUT1 "$ncol[3]\t$ncol[0]\t$ncol[1]\t$ncol[2]\n";
	    print OUT2 "$ncol[3]\t$ncol[0]\t$ncol[1]\t$ncol[2]\n";
	    $natom++;
	}else{
	    push @{$ch_coords[$n_chatom]}, @ncol;
	    $n_chatom++;
	}
}
continue{
    $i++; #print "\$i= $i\n";
}
#print "\$end\n";
print OUT1 "\$end\n\n";
print OUT2 "\$end\n\n";
### print external charge part
#print "\n\$external_charges\n";
print OUT2 "\n\$external_charges\n";
for($j=0;$j<$n_chatom;$j++){
#    print join("  ",@{$ch_coords[$j]}),"\n";
    for($k=0;$k<3;$k++){
	#print $ch_coords[$j][$k],"\t";
	print OUT2 $ch_coords[$j][$k],"\t";
    }
    #print $ext_charge,"\n";
    print OUT2 $ext_charge,"\n";
}

#print "\$end\n\n";
print OUT2 "\$end\n\n";

close(IN);
close(OUT1);
close(OUT2);

