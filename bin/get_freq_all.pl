#!/usr/bin/perl

$qcfile=$ARGV[0];

@fname=split(/\./, $qcfile);
$ffreq=$fname[0].".freq";
open(IN, "<$qcfile");
open(OUT, "> $ffreq");


$i=0;
while(<IN>){
    if($_ =~ /Frequency/ ){
	chomp($_);
	@line=split(/\s+/, $_);
	if($line[0] eq ""){ shift(@line)}
	    for($j=1;$j<=$#line;$j++){
		print OUT $line[$j]."\n";
	    }
	#print $_;
    }
}

close(IN);
close(OUT);
