#!/usr/bin/perl
# written by Joonho Park
# read  mof-co2 POSCAR 
# with direct coordinate ./a.pl poscar D
# with cartesian coordinate and extract co2 in .6co2.pos and .6co2.xyz; ./a.pl poscar (C)

if($#ARGV < 0){
    print "Error: it needs input file ";
    exit(1);
}
$fin=$ARGV[0];
$ext_mol=$ARGV[1];
$new_suffix="extmol";

$fout=$fin.".$new_suffix";

# GET CO2 coordinates
open(IN,"< $fin");
open(OUT,"> $fout");

if($ext_mol=~/ /){
    @ext_mol=split(/ +/,$ext_mol);
}else{
    print "Error: Usage: Mol Type should be separate with blank quoted by \" and \" \n";
    exit(2);
}

while(@ext_mol){$ext_kind[$i]=shift(@ext_mol); $ext_num[$i]=shift(@ext_mol); $i++;}
$ext_kind=join("",@ext_kind);
$new_line_kind=join("  ",@ext_kind);
$new_line_num=join("  ",@ext_num);
$nt_ext=0;
for($i=0;$i<=$#ext_kind;$i++){
    print "Atom: $ext_kind[$i] Number: $ext_num[$i] \n";
    $nt_ext+=$ext_num[$i];
}

$test="OFF";	# atom copy test
$tag="OFF"; 	# atom copy tag
$ikind=0;

@out_line=();
$ntotal=0;	# number of total atoms
$i=0;		# line index
$iatom=0;	# atom index
@ntotal_coord=();
@natom_kind=();
### read three atoms's coordinates
LINE: while($line=<IN>){
    if($i<5){ next LINE;}
    chomp($line);
    @line=split(/\s+/,$line);
#    print $line."\n";
    if($line[0] eq ""){shift(@line);}
    if($i==5){
	$string=join("",@line);
        die "There is numeric in atom line $!\n" if($string =~ /[^a-zA-Z]/); 
        for($j=0;$j<=$#line;$j++){push(@natom_kind,$line[$j]);}
	$new_line[$i]=$new_line_kind."\n";
        next LINE;
    }
    if($i==6){
        for($j=0;$j<=$#line;$j++){
            $ntotal+=$line[$j];
            $natom_num[$j]=$line[$j];	# number of atoms of each kind
        }
	$new_line[$i]=$new_line_num."\n";
        next LINE;
    }
    next LINE if($line[0] =~ /^[sS]/);	# for selected dynamics
    
    if($i==7){
	if($line[0] =~ /^[cC]/) {
	    print "Error::The input file is supposed to be Direct coordinates\n"; exit 11;
	}elsif($line[0] =~ /^[dD]/){
	    $new_line[$i]=$line."\n";
	    print "Direct coordinate is OK.\n"; next LINE;
	
	}else{ print "format error in cC|dD with $line[0]\n"; exit 12;}
    }
    # ntotal is calculated before 
    if($i>=8){
	$ntotal_coord[$iatom]=$line;
	$iatom++;
	if($iatom<$ntotal) {next LINE;}
	else {last LINE;}
    }

} continue {
    if($i==0){  $new_line[$i]=join("  ",@extr_atom);}
    if($i<5) {  $new_line[$i]=$line;} 
    $i++;
}

$total_line=$i;
if($iatom == $ntotal){
    	print "iatom = $iatom from coordinate, ntotal = $ntotal from atom number \n";
	print "input file is closed\n";
}else{	print "atom line $iatom is not matched with atom number $ntotal\n"; exit 9;}

close(IN);

#    print $itest, $ikind, $natom_kind[$ikind], $tag,"\n";

for($i=0; $i < 8; $i++){
    print $new_line[$i];
    print OUT $new_line[$i];
}

$ipivot=0;

for($i=0,$j=0,$k_line=0; $i<=$#natom_kind; $i++){
    $ipivot+=$natom_num[$i];
    ### atoms are different, then skip
    if($natom_kind[$i] eq $ext_kind[$j]){
	$jpivot=$ipivot-$ext_num[$j];
        print "here $ipivot $jpivot $ext_num[$j]\n";
        for($k=$jpivot;$k<$ext_num[$j]+$jpivot;$k++){
	    print $ntotal_coord[$k]."\n";
	    print OUT $ntotal_coord[$k]."\n";
	}
	$j++;
    }
}
	
close(OUT);
