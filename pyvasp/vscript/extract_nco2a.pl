#!/usr/bin/perl
# written by Joonho Park
# read  mof-co2 POSCAR 
# with direct coordinate ./a.pl poscar D
# with cartesian coordinate and extract co2 in .6co2.pos and .6co2.xyz; ./a.pl poscar (C)

if($#ARGV < 0){
    print "Error: it needs input file ";
    exit(0);
}
$fin=$ARGV[0];

$nmol=$ARGV[1];

if($#ARGV >= 1){
    $format = $ARGV[2];
}else{
    $format = "C";
}

if($format =~ /^[dD]/){
    $format = "dirt";
}elsif($format =~ /^[cC]/){
    $format = "cart"
}

$nmol=6 if(!defined($nmol));

$new_suffix="${nmol}co2.$format";
$xyz_suffix="${nmol}co2.xyz";
if($fin =~ /\./){
    @name=split(/\./,$fin);
    $fout=$name[0].".$new_suffix";
    $fxyz=$name[0].".$xyz_suffix";
}else{
    $fout=$fin.".$new_suffix";
    $fxyz=$fin.".$xyz_suffix";
}

if($format =~ /^[dD]/){
    print "Extract 6CO2 in direct coordinate in $fout\n";
    $xyztag="OFF";
}elsif($format =~ /^[cC]/){
    print "Extract 6CO2 in cartesian coordinate in $fout and xyz format in $fxyz\n";
    open(XYZ,"> $fxyz");
    $xyztag="ON";
}else{ print "Error:: what's the format or coordination\n"; exit;}

# GET CO2 coordinates
open(IN,"< $fin");
open(OUT,"> $fout");

@extr_atom=qw( O C );
@natom_ea=qw( 2 1 );

for($i=0;$i<=$#extr_atom;$i++){
    $natom_extracted[$i]=$nmol*$natom_ea[$i];
}

$nt_extr=$nmol*3;
$test="OFF";	# atom copy test
$tag="OFF"; 	# atom copy tag
$ikind=0;

@out_line=();
$ntatoms=0;
$i=0;
### read three atoms's coordinates
LINE: while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($i==0){  print OUT join("  ",@extr_atom), "\n";
		push(@out_line,join("  ",@extr_atom)."\n"); next LINE;}
    if($i==1){ print OUT $line."\n"; 
		push(@out_line,$line."\n"); next LINE;}
    if($i<5){  print OUT $line."\n"; 
		push(@out_line,$line."\n"); next LINE;}
    if($i==5){
        if($line[0] !~ /[A-Z]/) {
            print "Error in $fin format: there are not atom symbols in $i+1 line\n";
            exit(0);
        }else{ for($j=0;$j<=$#line;$j++){push(@atoms,$line[$j]);} }
	print OUT join("  ",@extr_atom)."\n";
	push(@out_line,join("  ",@extr_atom)."\n");
	if($xyztag eq "ON"){
	    print XYZ $nt_extr,"\n";
	    print XYZ join("  ",@extr_atom)."\n";
	}
        next LINE;
    }
    if($i==6){
        for($j=0;$j<=$#line;$j++){
            $ntatoms+=$line[$j];
            $natom_kind[$j]=$line[$j];	# number of atoms of each kind
        }
	print OUT join("  ",@natom_extracted),"\n";
	push(@out_line,join("  ",@natom_extracted)."\n");
        next LINE;
    }
    next LINE if($line[0] =~ /^[sS]/);	# for selected dynamics
    if($line[0] =~ /^[cC]/) {
	if($xyztag eq "ON"){
	    print OUT $line."\n"; next LINE;
	}else{print "Error::The input file is supposed to be Cartesian coordinates\n"; exit;}
    }elsif($line[0] =~ /^[dD]/){
	if($xyztag eq "OFF"){
	    print OUT $line."\n"; next LINE;
	}else{print "Error::Input file has direct coordinates\n"; exit;}
    }

    # AFTER knowing @atoms, $natom_kind, we calculate the pivot
    # WRITE or SKIP test 
    if($test eq "OFF" and $itest ==0){
    	EXTR_LIST: for($j=0; $j<=$#extr_atom; $j++){
	    if($atoms[$ikind] eq $extr_atom[$j]){	# start from $ikind=0, first kind
		$tag="ON"; $iextr=$j;			# $iextr saves the kind to be extracted
		last EXTR_LIST;
	    }
	}
	$test = "ON";
	$itest = 0; 	# start number of kind-@atoms
    }

    # Test for atom check
#    if($tag_test eq "OFF"){
#	for($j=0;$j<=$#atoms;$j++){
#	    if($extr_atom[$j] eq $atom[0]){
#		$tag_extr="ON";
#		last;
#	    }else{
#		$tag_extr="OFF";
#	    }
#	}
#	if($tag_extr eq "ON"){
#	    $n_extr=
#	$tag_test="ON";
#	    

#    print $itest, $ikind, $natom_kind[$ikind], $tag,"\n";
    if($itest < $natom_kind[$ikind]){
	if($tag eq "ON" and $natom_kind[$ikind]-$natom_extracted[$iextr] <= $itest ){
	    print "OUT: $ikind  $natom_kind[$ikind] $natom_extracted[$iextr] $itest \n";
	    print OUT "$line[0]  $line[1]  $line[2]\n";
	    print XYZ "$extr_atom[$iextr]  $line[0]  $line[1]  $line[2]\n" if($xyztag eq "ON");
	}
	$itest++;
	if($itest == $natom_kind[$ikind]){
	    $ikind++; $test  = "OFF"; $tag   = "OFF"; $itest = 0;
	} 
    }
} continue {
    $i++;
}

open(OUTTest,">out_test");
for($i=0;$i<=$#out_line;$i++){
    print OUTTest $out_line[$i];
}

close(OUTTest);


close(IN); close(OUT); close(XYZ);
