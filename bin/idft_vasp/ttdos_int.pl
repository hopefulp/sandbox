#!/usr/bin/perl

$fin=$ARGV[0];

open(IN,"<$fin1");

@ene=();

for($f=0;$f<2;$f++){
    open(IN,"<@ARGV[$f]");
    $i=0;
    while($line=<IN>){
        chomp($line);
        @line=split(/\s+/,$line);
        if($line[0] eq ""){shift(@line);}
 
        $ene[$f][$i]=$line[0];
        $dos[$f][$i]=$line[1];
    }continue{
    	$i++;
    }
    if($f==0){    $niene_1=$i;}
    if($f==1){    $niene_2=$i;}

    close(IN);
}

for($i=0;$i<$niene;$i++){
    print $ene[0][$i],"\t",$dos[0][$i],"\n";
}


