#!/usr/bin/perl
# written by Joonho Park
# read file, find key_word, print line

$fin=shift;
$nkey=shift || 1;
$nline=shift || 10;

open(IN,"<$fin");

$i=0;
$iline=0;
$ikey=0;
$tag="OFF";

LINE: while($line=<IN>){
#    $i++;
#    print $line;

    if($line =~ /FREQ\]/ and $tag eq "OFF"){
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
    

