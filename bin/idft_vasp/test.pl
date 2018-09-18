#!/usr/bin/perl

# written by Joonho Park
# read POSCAR or CONTCAR format and write xyz file


if($#ARGV < 0){
    print "Error: input POSCAR's ";
    exit(0);
}

$fin=shift(@ARGV);
if($fin =~ /\./){
    @name=split(/\./,$fin);
    $fout=$name[0].".xyz";
}else{
    $fout=$fin.".xyz";
}
@weird_coord=();
for($i=0;$#ARGV >= 0;$i++){
    $weird_coord[$i]=shift(@ARGV);
}

print "OUTFILE:: $fout\n";

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

    next LINE if($i==0);
    if($i==1){ $scale=$line; next LINE;}
    if($i<5){
	for($j=0;$j<3;$j++){$lattice_vector[$i-2][$j]=$line[$j];}
	next LINE;
    }
    if($i==5){
	if($line[0] !~ /[A-Z]/) {
	    print "Error in $fin format: there are not atom symbols in 6th line\n";
	    exit(0);
	}else{ for($j=0;$j<=$#line;$j++){push(@atoms,$line[$j]);} }
	next LINE;
    }
    if($i==6){
	for($j=0;$j<=$#line;$j++){ 
	    $ntatoms+=$line[$j];
	    $nsum[$j]=$ntatoms;
 	}
	next LINE;
    }
    next LINE if($line[0] =~ /[a-zA-Z]/);
    
    for($j=0;$j<3;$j++){ $coord[$iatom][$j]=$line[$j];}
    $iatom++;
    if($iatom >= $ntatoms){ last;}		
} continue {
    $i++;
}

close(IN); 

print $ntatoms."\n";
print OUT $ntatoms."\n";
print join(" ",@atoms)."\n";
print OUT join(" ",@atoms)."\n";

# usage: a.pl CONTCAR 30 34
# Cr 30(O) imov=2(z) -1; 34(O) imov=0 -1: @imov=(2,0);  @jmov=(-1,-1);
# Cu 28 30 34				: @imov=(2,2,0);@jmov=(1,-1,-1);
# Fe 30					: @imov=(2); 	@jmov=(-1);
@imov=(2);
@jmov=(-1);
for($i=0;$i<$ntatoms;$i++){
    if($i<$nsum[0]){
#	print "$atoms[0]  ";
	print OUT "$atoms[0]  ";
    }elsif($i<$nsum[1]){
#	print "$atoms[1]  ";
	print OUT "$atoms[1]  ";
    }elsif($i<$nsum[2]){
#	print "$atoms[2]  ";
	print OUT "$atoms[2]  ";
    }else{
#	print "$atoms[3]  ";
	print OUT "$atoms[3]  ";
    }

    if($#weird_coord >= 0){
	for($j=0;$j<=$#weird_coord;$j++){
	    if($i==$weird_coord[$j]-1){
	        $imov=shift(@imov);
		$jmov=shift(@jmov);
		$orig=$coord[$i][$imov];
		$coord[$i][$imov]=$orig+$jmov;
		
		print "orig $orig changed $coord[$i][$imov]\n";
		shift(@weird_coord);
		last;
	    }
	}
    }

    ($x,$y,$z)=&xyz($coord[$i][0],$coord[$i][1],$coord[$i][2]);
    
#    print $x*$scale,"\t",$y*$scale,"\t",$z*$scale,"\n";
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
