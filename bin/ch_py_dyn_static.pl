#!/usr/bin/perl
# written by Joonho Park
# read file, find key_word, print line

$fin=shift @ARGV;

open(IN,"<$fin");
@fname=split(/\./, $fin);
$fname=$fname[0]."_static.py";
open(OUT,">$fname ");

#print $fname;

$i=0;
$iline=0;
$ikey=0;
$tag="OFF";
$key1="while";

LINE: while($line=<IN>){
#    $i++;
#    print $line;

    if($line =~ /$key1/){
        
	    $ikey++;
	print $line;
	if($tag eq "OFF"){
	    print "found key: $ikey times\n";
	    if($ikey == $nkey){	    $tag="ON";}
	}
	next LINE;
    }elsif($tag eq "ON"){
	print $line;
	$iline++;
	if($iline == $nline){
	    $tag="OFF";
	}
	next LINE;
    }
}
continue{
    $i++;
}
close(IN);
    

