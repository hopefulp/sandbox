#!/usr/bin/perl
# written by Joonho Park
# if read a.mol, write a.xyz
# if read a.xyz, write a.mol 

$fin=shift @ARGV;
$f_charge=shift @ARGV;
@list_atom=@ARGV;

$charge=0;
$multi=1;

#use lib '/qcfs/joonho/modules';
use lib '/home/joonho/sandbox_gl/mod_perl';
use IOFile2;

# output file
($fout,$suffix)=IOFile2::file2($fin);

if(0<=$#ARGV){ $charge=shift @ARGV;}
if(0<=$#ARGV){ $multi=shift @ARGV;}

#print "charge $charge multi $multi\n";

open(IN,"<$fin");
open(OUT,">$fout");

@f_line=();
### read three atoms's coordinates
$i=0;
$natom=0;
while($line=<IN>){
    # skip three lines 
    if($i < 2){next;}
    else{	
	@ncol=split(/\s+/,$line);
	if($ncol[0] eq "") {shift @ncol;}
	#print join("\t",@ncol), "\n";
	# if infile is mol-file, don't copy last line of "$end"
	# number of copied lines are the number of atoms
	if( $#ncol == 3 ){
	    push(@f_line,$line);
	    $natom++;
	}
    }
}
continue{
    $i++;
}

# 2 lines pre-script for outfile 
if($suffix eq "mol"){
    print OUT "\$molecule\n";
    print OUT "\t$charge\t$multi\n";
    print     "\$molecule\n";
    print     "\t$charge\t$multi\n";
}elsif($suffix eq "xyz"){
    print OUT $natom."\n";
    print OUT $fout."\n";
    print     $natom."\n";
    print     $fout."\n";
}
# write atom list
print OUT @f_line;
print     @f_line;

# closing line for mol-file
if($suffix eq "mol"){
    print OUT "\$end\n";
    print     "\$end\n";
}

close(IN); 
close(OUT);
