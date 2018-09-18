#!/usr/bin/perl
# written by Joonho Park

$dir=$ARGV[0];

open(IN,"<$dir/DOSCAR");

$i=0;
$iatom=1;
$sum=0;
LINE: while($line=<IN>){
#    print $line;
    next if($i<5);
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($i==5){
        $delimit_atom=$line;
        $Emax=$line[0];
        $Emin=$line[1];
        $Nene=$line[2];
        $ene_Fermi=$line[3];
        $dnknow=$line[4];
        next LINE;
    }
  
} continue {
    $i++;
}
close(IN);

if($#ARGV == 1){
    $dir=$ARGV[1];
    open(IN,"<$dir/DOSCAR");

    $i=0;
    $iatom=1;
    $sum=0;
    LINE2: while($line=<IN>){
#    print $line;
        next if($i<5);
        chomp($line);
        @line=split(/\s+/,$line);
        if($line[0] eq ""){shift(@line);}
 
        if($i==5){
            $delimit_atom=$line;
            $Emax2=$line[0];
            $Emin2=$line[1];
            $Nene2=$line[2];
            $ene_Fermi2=$line[3];
            $dnknow2=$line[4];
            next LINE2;
        }

    } continue {
    	$i++;
    }
    close(IN);
}



print $ene_Fermi,"\t",$ene_Fermi2 ,"\t", $ene_Fermi2-$ene_Fermi,"\n";


