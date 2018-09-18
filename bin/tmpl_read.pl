#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

@pinch=(1, 25, 26, 51);

$fin=$ARGV[0];

@fname=split(/\./,$fin);

$sys=$fname[0];
$ext=$fname[1];
if($fname[1] ne "xyz"){
    print "input error: the suffix should be xyz of qchem input file\n";
    exit(0);
}

$fout="$sys.mol";

open(IN,"<$fin");
open(OUT,">$fout");

$line1="\$molecule";
$line2="\$end";
@f_line=();
### read three atoms's coordinates
$i=0;
$j=0;
$p=$pinch[$j];

while($line=<IN>){
    if($i==$p+1){
	push(@f_line,$line);
    	@field=split(/\s+/,$line);
    	if($field[0] eq "") {shift(@field);}
	$p=$pinch[++$j];
    }else{

    }
}continue{
    $i++;
}

print $line1,"\n";
print @f_line;
print $line2,"\n";
print OUT $#f_line+1, "\n\n", @f_line;

close(IN); 
close(OUT);
