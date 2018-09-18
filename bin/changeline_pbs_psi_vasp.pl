#!/usr/bin/perl
# written by Joonho Park
# read a.csh and write t.csh

use Switch;

if($#ARGV < 0){
    print "Usage:: $0 pbs-psi-cell.csh dirname vasp_exe\n";
    exit;
}elsif($#ARGV == 1){
    $fin="pbs-psi-cell.csh";
    $dir=$ARGV[0];
    $exe_vasp=$ARGV[1];
}elsif($#ARGV == 2){
    $fin=$ARGV[0];
    $dir=$ARGV[1];
    $exe_vasp=$ARGV[2];
}

$exe_string="export EXEC=\"/opt/applic/vasp/bin/vasp-5.3.3-xe11-static-openmpi-1.6.3-";
$exe_string_xy="export EXEC=\"/opt/applic/vasp/bin/xy-vasp-5.3.3-xe11-static-openmpi-1.6.3-";
switch ($exe_vasp){
    case "gamma" { $new_string=$exe_string."gamma\"\n"; }
    case "full"  { $new_string=$exe_string."full\"\n"; }
    case "half"  { $new_string=$exe_string_xy."half\"\n"; }
    else 	 { print "Usage:: there is no vasp execution as 2nd arg $exe_vasp"; exit; }
}

$fout="t.sh";
print "Input file is $fin\n";

open(IN,"< $fin") || die "There is not pbs script as input file \"$fin\"\n";
open(OUT,"> $fout");
#close(IN); close(OUT);
#exit;

@f_line=();
### read three atoms's coordinates
$i=0;
LINE: while($line=<IN>){
#    $i++;
    @line=split(/\s+/,$line);
    if($line[0] eq "") { shift @line; }
#    print "$line[0]\t\t$line[1]\n";
    if($line[1] eq "-N"){
	print OUT "#PBS -N $dir\n";
    }elsif($line[0] =~ /^EXEC/){
	print OUT $new_string;
	print $new_string;
#	print OUT "set EXEC = \"/opt/applic/vasp/bin/vasp-5.3.3-xe11-static-openmpi-1.6.3-$exe_vasp\"\n";	
    }else{
	print OUT $line;
    }
}
continue{
    $i++;
}
close(IN); 
close(OUT);

system("qsub $fout");
print "The job will be running on 1 nodes with 16 processes\n";

