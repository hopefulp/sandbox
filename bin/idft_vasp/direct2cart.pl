#!/usr/bin/perl

# written by Joonho Park
# read POSCAR or CONTCAR format in fractional coordinate and write in cartesian coordinates

if($#ARGV < 0){
    print "Error: input POSCAR's ";
    exit(0);
}

$fin=$ARGV[0];
$new_suffix="cart";
if($fin =~ /\./){
    @name=split(/\./,$fin);
    $fout=$name[0].".$new_suffix";
}else{
    $fout=$fin.".$new_suffix";
}

open(IN, "<$fin");
open(OUT, ">$fout");
$i=0;
@atoms=();
@nsum=();
$ntatoms=0;
$iatom=0;
LINE: while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($i==0){ print OUT $line."\n"; next LINE;}
    if($i==1){ $scale=$line; print OUT "    1.0\n"; next LINE;}
    if($i<5){
	for($j=0;$j<3;$j++){
	    $lattice_vector[$i-2][$j]=$line[$j]*$scale;
	    print OUT "  ",$lattice_vector[$i-2][$j];
	}
	print OUT "\n";
	next LINE;
    }
    if($i==5){
	if($line[0] !~ /[A-Z]/) {
	    print "Error in $fin format: there are not atom symbols in 6th line\n";
	    exit(0);
	}else{ for($j=0;$j<=$#line;$j++){push(@atoms,$line[$j]);} print OUT join(" ",@atoms)."\n"; }
	next LINE;
    }
    if($i==6){
	for($j=0;$j<=$#line;$j++){
	    print OUT "  $line[$j]";
	    $ntatoms+=$line[$j];
	    $nsum[$j]=$ntatoms;
 	}
	print OUT "\n";
	next LINE;
    }
    if($line[0] =~ /^[sS]/) {print OUT $line."\n";    next;}
    if($line[0] =~ /^[dD]/) {print OUT "Cartesian\n"; next;}
    for($j=0;$j<3;$j++){ $coord[$iatom][$j]=$line[$j];}
    $iatom++;
    if($iatom >= $ntatoms){ last;}		
} continue{
    $i++;
}
close(IN); 

print "total number of atom :: $iatom\n";


for($i=0;$i<$ntatoms;$i++){
    ($x,$y,$z)=&xyz($coord[$i][0],$coord[$i][1],$coord[$i][2]);
    print $x*$scale,"\t",$y*$scale,"\t",$z*$scale,"\n";
    print OUT $x*$scale,"\t",$y*$scale,"\t",$z*$scale,"\n";
#    exit;
}

sub xyz{
    my ($ap,$bp,$cp)=@_;
    my ($x, $y, $z);
    $x=$lattice_vector[0][0]*$ap+$lattice_vector[1][0]*$bp+$lattice_vector[2][0]*$cp;
    $y=$lattice_vector[0][1]*$ap+$lattice_vector[1][1]*$bp+$lattice_vector[2][1]*$cp;
    $z=$lattice_vector[0][2]*$ap+$lattice_vector[1][2]*$bp+$lattice_vector[2][2]*$cp;

    return ($x,$y,$z);
}
close(OUT);
