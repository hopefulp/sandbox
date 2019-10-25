#!/usr/bin/perl
# written by Joonho Park
# read  mof-co2 POSCAR and extract co2 in .6co2.xyz

if($#ARGV < 0){
    print "Error: input POSCAR with Cartesian coordinates ";
    exit(0);
}

$fin=$ARGV[0];
$new_suffix="6co2.pos";
$xyz_suffix="6co2.xyz";
if($fin =~ /\./){
    @name=split(/\./,$fin);
    $fout=$name[0].".$new_suffix";
    $fxyz=$name[0].".$xyz_suffix";
}else{
    $fout=$fin.".$new_suffix";
    $fxyz=$fin.".$xyz_suffix";
}

print "New file is $fout; xyz format is in $fxyz\n";

# GET CO2 coordinates
open(IN,"< $fin");
open(OUT,"> $fout");
open(XYZ,"> $fxyz");

@extr_atom=qw( O C );
@extr_natom=qw( 12 6 );

$nt_extr=18;
$xyztag="ON";
$test="OFF";	# atom copy test
$tag="OFF"; 	# atom copy tag
$ikind=0;
### read three atoms's coordinates
LINE: while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($i==0){ print OUT join("  ",@extr_atom), "\n"; next LINE;}
    if($i==1){ print OUT $line."\n"; next LINE;}
    if($i<5){  print OUT $line."\n"; next LINE;}
    if($i==5){
        if($line[0] !~ /[A-Z]/) {
            print "Error in $fin format: there are not atom symbols in 6th line\n";
            exit(0);
        }else{ for($j=0;$j<=$#line;$j++){push(@atoms,$line[$j]);} }
	print OUT join("  ",@extr_atom)."\n";
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
	print OUT join("  ",@extr_natom),"\n";
        next LINE;
    }
    next LINE if($line[0] =~ /^[sS]/);
    if($line[0] =~ /^[cC]/) {print OUT $line."\n"; next LINE;}
    if($line[0] =~ /^[dD]/) {print "Error::Input file has direct coordinates\n"; exit;}
    # write or skip test 
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
#    print $itest, $ikind, $natom_kind[$ikind], $tag,"\n";
    if($itest < $natom_kind[$ikind]){
	if($tag eq "ON" and $natom_kind[$ikind]-$extr_natom[$iextr] <= $itest ){
	    print "OUT: $ikind  $natom_kind[$ikind] $extr_natom[$iextr] $itest \n";
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

close(IN); close(OUT); close(XYZ);
