#!/usr/bin/perl
# written by Joonho Park
# read a.csh and write t.csh

#use List::MoreUtils 'true';

$dir=$ARGV[0];
$fname="CHGCAR";
$fname1="CHGCAR1";
open(IN,"<$dir/$fname") || die "There is no $fname" ;
open(OUT,">$dir/$fname1");

### read three atoms's coordinates
$i=0;
$mod_tag="OFF";
LINE: while($line=<IN>){
#    $i++;
    if($i == 4 || $i == 5 ){
      	@line=split(/s+/, $line);
        if($line[0] eq "") { shift @line; }
    	if($line !~ /\d/){
	    $mod_tag="ON";
	    print "Atom symbol line in CHGCAR was deleted\n";
	    next LINE;
	}
#	print $i, ":", $line, @bar1;
    }
    print OUT $line;
}
continue{
    $i++;
}

close(IN); 
close(OUT);

### if CHGCAR is not normal file die
if($i < 10){
    print "CHGCAR has not proper file with $i lines\n";
    exit;
}

if($mod_tag eq "ON"){
    system("mv $dir/$fname1 $dir/$fname");
}else{
    print "Atom symbol line in CHGCAR was not detected\n";
}


