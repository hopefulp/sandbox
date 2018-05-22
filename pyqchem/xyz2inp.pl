#!/usr/bin/perl
# written by Joonho Park
# read a.xyz or stdin and write a.mol

$tag_infile="";
$charge=0;
$multi=1;

if(0 <= $#ARGV ){
    $xyz=$ARGV[0];
    @fname=split(/\./,$xyz);
    $tag_infile="T";
    if($ARGV[1] ne ""){
	$charge=$ARGV[1];
	$multi=$ARGV[2];
    }
}

if($tag_infile ne ""){
    $fout="$fname[0].mol";
    open(IN,"<$xyz");
}else{
#    IN="STDIN";
    $fout="t.mol";
}
open(OUT,">$fout");


### read three atoms's coordinates
$i=0;

print OUT "\$molecule\n\t$charge\t$multi\n";

if($tag_infile ne ""){
    while($line=<IN>){
        $i++;
        if($i <=2){ next;}
        #print $line;
        print OUT $line;
    } 
}else{
    while($line=<STDIN>){
        $i++;
        if($i <=2){ next;}
        print $line;
        print OUT $line;
    } 
}

print OUT "\$end\n";

close(IN); 
close(OUT);
