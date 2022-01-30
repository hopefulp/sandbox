#!/usr/bin/perl
# written by Joonho Park
# read file, find key_word, print line

$fin=shift;
$nkey=shift || 1;
$nline=shift || 10;

if($fin eq ""){
    print "Usage:: $0 qchem_outfile | awk '{print \$3}'\n";
    exit(0);
}    

open(IN,"<$fin");

$key_word="Valence    Rydberg";
$end_word="=====";

$i=0;
$iline=0;
$ikey=0;
$tag="OFF";


LINE: while($line=<IN>){
#    $i++;
#    print $line;

    if($line =~ /$key_word/ and $tag eq "OFF"){
	    $ikey++;
	    print $line;
	    if($tag eq "OFF"){
	        print "found key: $ikey times\n";
            ### print 3 times, total, alpha, beta
            $tag="ON";
	        #if($ikey == $nkey){	    $tag="ON";}
	    }
	    next LINE;
    }elsif($tag eq "ON"){
	    #if($iline == $nline){
        if($line =~ /$end_word/){
	        $tag="OFF";
            next LINE;
	    }   
	    print $line;
	    $iline++;
	    next LINE;
    }
}
continue{
    $i++;
}
close(IN);
    

