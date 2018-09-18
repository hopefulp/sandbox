#!/usr/bin/perl
# written by Joonho Park
# read INCAR for SCF and write INCAR for DOS

$fin=$ARGV[0];
$fout="INCAR.dos";

if($#ARGV < 0){
    print "Usage:: $0 INCAR \n";
    print "   make $fout\n";
    exit(1);
}

open(IN,"<$fin");
open(OUT,">$fout");

$tag_algo="N";
$tag_nsw="N";
$tag_npar="N";
$tag_ib="N";
$i=0;
while($line=<IN>){
    @line=split(/\W+/,$line);
    if($line[0] eq "") {shift @line;}
#    print $line[0],"\n";

    if($line =~ /ISTART/){	print OUT  "ISTART = 0\n"; 	next; }
    if($line =~ /ICHARG/){	print OUT  "ICHARG = 11\n"; 	next; }
    if($line[0] eq "NSW" and $tag_nsw eq "N")  	{print OUT  "NSW = 0\n"; $tag_nsw="Y";next; }
    if($line =~ /IBRION/ and $tag_ib eq "N"){	print OUT  "IBRION = -1\n"; $tag_ib="Y"; next; }
    if($line[0] eq "ALGO" and $tag_algo eq "N") {print OUT  "ALGO = Normal\n"; $tag_algo="Y"; next; }
    if($line =~ /NEDOS/) {	print OUT  "NEDOS = 4000\n";	 next; }
    if($line =~ /NPAR/ and $tag_npar eq "N")	{print OUT  "NPAR = 1\n"; $tag_npar="Y"; next;}
    if($line =~ /LAECHG/){	print OUT  "LAECHG = .FALSE.\n"; next;}
    if($line =~ /LWAVE/) {	print OUT  "LWAVE  = .TRUE.\n";  next;}
    if($line =~ /LCHARG/){	print OUT  "LCHARG = .FALSE.\n"; next;}
    print OUT  $line;
}
continue{
    $i++;
}
close(IN); 
close(OUT);


