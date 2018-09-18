#!/usr/bin/perl
# written by Joonho Park
# read INCAR for SCF and write INCAR for DOS

$fin=$ARGV[0];
$fout="INCAR.pchg";

if($#ARGV < 0){
    print "Usage:: $0 INCAR \n";
    print "   make $fout\n";
    exit(1);
}

open(IN,"<$fin");
open(OUT,">$fout");

$tag_algo="N";
$tag_nsw="N";
#$tag_npar="N";		# remove NPAR
$tag_ib="N";
$tag_potim="N";

### for vibration calculation
$tag_write="N";
$tag_nfree="N";

### for pchg calculation
$tag_iband="N";

$i=0;
while($line=<IN>){
    @line=split(/\W+/,$line);
    if($line[0] eq "") {shift @line;}
#    print $line[0],"\n";

    if($line =~ /ISTART/){	print OUT  "ISTART = 1\n"; 	next; }
    if($line =~ /ICHARG/){	print OUT  "ICHARG = 11 # USE non-self-consistent cal.\n"; 	next; }
    if($line[0] eq "NSW" and $tag_nsw eq "N")  	{print OUT  "#NSW = 0 \n"; $tag_nsw="Y";next; }
    if($line =~ /IBRION/ and $tag_ib eq "N"){	print OUT  "#IBRION = 2\n"; $tag_ib="Y"; next; }
    if($line[0] eq "ALGO" and $tag_algo eq "N") {print OUT  "ALGO = Fast\n"; $tag_algo="Y"; next; }
    if($line =~ /NPAR/ and $tag_npar eq "N")	{print OUT  "NPAR = 1\n"; $tag_npar="Y"; next;}
    if($line =~ /POTIM/ and $tag_potim eq "N") { print OUT  "#POTIM = 0.01\n"; $tag_potim="Y"; next;}
    if($line =~ /NWRITE/ and $tag_write eq "N") {print OUT  "#NWRITE = 3\n"; $tag_write="Y"; next;}
    if($line =~ /LAECHG/){	print OUT  "LAECHG = .FALSE.\n"; next;}
    if($line =~ /LWAVE/) {	print OUT  "LWAVE  = .FALSE.\n";  next;}
    if($line =~ /LCHARG/){	print OUT  "LCHARG = .FALSE.\n"; next;}
    if($line =~ /IBAND/) {	print OUT  "IBAND = 1\n"; $tag_iband="Y"; next;}
    if($line =~ /KPUSE/) {	print OUT "KPUSE = 1 2 3 4 5 6 7 8 \n"; next;}
    if($line =~ /LSEPB/) {	print OUT "LSEPB = .TRUE.\n"; next;}
    if($line =~ /LSEPK/) {	print OUT "LSEPK = .TRUE.\n"; next;}
    if($line =~ /LPARD/) {	print OUT "LPARD = .TRUE.\n"; next;}

### To remove DOS
    if($line =~ /NEDOS/) {	print OUT "#".$line; next;}

    print OUT  $line;
}
continue{
    $i++;
}

if($tag_iband eq "N"){
      print OUT  "IBAND = 1\n";
      print OUT "KPUSE = 1\n"; 
      print OUT "LSEPB = .TRUE.\n"; 
      print OUT "LSEPK = .TRUE.\n"; 
      print OUT "LPARD = .TRUE.\n"; 
}

close(IN); 
close(OUT);


