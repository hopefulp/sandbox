#!/usr/bin/perl
# written by Joonho Park
# read Me.xyz and write Me_co2.mol

@pinch=qw( 1 2 3 4 5 6 25 26 27 28 29 30 31 32 33 34 35 36 61 62 63 64 65 66);

for($i=0,$k=0;$i<6;$i++){
    for($j=0;$j<4;$j++){
	$pair[$i][$j]=


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

@atom_list=qw(Mg Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn);

$fin=$ARGV[0];
@fname=split(/\./,$fin);
$Me=$fname[0];
$ext=$fname[1];
if($fname[1] ne "xyz"){
    print "input error: the suffix should be xyz of qchem input file\n";
    exit(0);
}

# Find metal location
for($i=0;$i<=$#atom_list;$i++){
    if($atom_list[$i] eq $Me){
	$iMe=$i;
	last;
    }
}
#print "iMe = $iMe\n";
$fcharge="metal_pcharge.dat";
open(IN,"<$fcharge");
$i=0;
while(<IN>){
    @me_charges=split(/\s+/,$_);
    if($me_charges[0] eq ""){
	shift(@me_charges);
    }
    for($j=0;$j<=$#me_charges;$j++){
	$all_charges[$i][$j]=$me_charges[$j];
    }
#    print "Me charges: $me_charges[$j-1], $all_charges[$i][$j-1]\n";
}continue{
    $i++;
}
close(IN);

$fout="${Me}_co2_image.mol";
$fout2="${Me}_co2.mol";
$chg=$chg{$Me};
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
	if($p<=6){ push(@ext_chg,$line);}
	else	 { push(@mol_line,$line);}
#    	@field=split(/\s+/,$line);
#    	if($field[0] eq "") {shift(@field);}
	$p=$pinch[++$j];
    }else{

    }
}continue{
    $i++;
}

$all_ext=join("  ",@ext_chg);
@field=split(/\s+/,$all_ext);
if($field[0] eq "") {shift(@field);}

#print @field,"\n";

print $line1,"\n";
print "  0  1\n";
print @mol_line;
print $line2,"\n";
print "\$external_charges\n";
for($i=0,$j=0;$i<=$#field;$i+=4,$j++){
    {	print "  $field[$i+1]  $field[$i+2]  $field[$i+3] $all_charges[$iMe][$j]\n"; }
}
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
for($i=0,$j=0;$i<=$#field;$i+=4,$j++){
    {	print OUT"  $field[$i+1]  $field[$i+2]  $field[$i+3] $all_charges[$iMe][$j]\n"; }
}
print OUT "\$end\n";

close(IN); 
close(OUT);close(OUT2);
