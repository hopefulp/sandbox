#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

@pinch=(1, 25, 26, 51);
%chg=( 
'Mg'=>1.683 ,
'Ca'=>1.610 ,
'Sc'=>1.995 ,
'Ti'=>1.659 ,
'V'=>1.589 ,
'Cr'=>1.398 ,
'Mn'=>1.467 ,
'Fe'=>1.354 ,
'Co'=>1.274 ,
'Ni'=>1.215 ,
'Cu'=>1.077 ,
'Zn'=>1.368 
);


$fin=$ARGV[0];

@fname=split(/\./,$fin);

$sys=$fname[0];
$ext=$fname[1];
if($fname[1] ne "xyz"){
    print "input error: the suffix should be xyz of qchem input file\n";
    exit(0);
}

$fout="${sys}_co2.mol";
$fout2="${sys}__co2.mol";
$chg=$chg{$sys};
#print $chg,"\n";
open(IN,"<$fin");
open(OUT,">$fout");
open(OUT2,">$fout2");
$line1="\$molecule";
$line2="\$end";
@mol_line=();
@ext_chg=();
### read three atoms's coordinates
$i=0;
$j=0;
$p=$pinch[$j];

while($line=<IN>){
    if($i==$p+1){
	if($p==1){ push(@ext_chg,$line);}
	else	 { push(@mol_line,$line);}
#    	@field=split(/\s+/,$line);
#    	if($field[0] eq "") {shift(@field);}
	$p=$pinch[++$j];
    }else{

    }
}continue{
    $i++;
}

@field=split(/\s+/,$ext_chg[0]);
if($field[0] eq "") {shift(@field);}

print $line1,"\n";
print "  0  1\n";
print @mol_line;
print $line2,"\n";
print "\$external_charges\n";
print " $field[1]  $field[2]  $field[3] $chg\n";
print "\$end\n";

print OUT $line1,"\n";
print OUT2 $line1,"\n";
print OUT "  0  1\n";
print OUT2 "  0  1\n";
print OUT @mol_line;
print OUT2 @mol_line;
print OUT $line2,"\n";
print OUT2 $line2,"\n";
print OUT "\$external_charges\n";
print OUT " $field[1]  $field[2]  $field[3] $chg\n";
print OUT "\$end\n";

close(IN); 
close(OUT);close(OUT2);
