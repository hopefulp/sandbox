#!/usr/bin/perl

# written by Joonho Park
# read POSCAR or CONTCAR format in fractional coordinate and write in cartesian coordinates

if($#ARGV < 0){
    print "Error: input POSCAR's ";
    exit(0);
}

$fin=$ARGV[0];
#$new_suffix="cart";
$new_suffix="ct";
if($fin =~ /\./){
    @name=split(/\./,$fin);
    $fout=$name[0].".$new_suffix";
}else{
    $fout=$fin.".$new_suffix";
}

open(IN, "<$fin");
open(OUT, ">$fout");
print "write to $fout\n";
$i=0;
@atoms  = ();
@nsum   = ();
$ntatoms= 0;
$iatom  = 0;
$tag_sel    = 0;
LINE: while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($i==0){ 
        print $line."\n"; 
        print OUT $line."\n"; 
        next LINE;}
    if($i==1){
        $scale=$line; 
        print OUT "    1.0\n"; 
        print "    1.0\n"; 
        next LINE;}
    if($i<5){
	    for($j=0;$j<3;$j++){
            $lattice_vector[$i-2][$j]=$line[$j]*$scale;     # $i-2 to make it 0
            printf OUT  "\t%10.5f", $lattice_vector[$i-2][$j];
            printf      "\t%10.5f", $lattice_vector[$i-2][$j];
	    }
        print OUT   "\n";
        print       "\n";
        next LINE;
    }
    if($i==5){
        if($line[0] !~ /[A-Z]/) {
            print "Error in $fin format: there are not atom symbols in 6th line\n";
            exit(0);
        }else{ 
            for($j=0;$j<=$#line;$j++){
                push(@atoms,$line[$j]);} 
            print OUT   join(" ",@atoms)."\n"; 
            print       join(" ",@atoms)."\n"; 
        }
        next LINE;
    }
    if($i==6){
	for($j=0;$j<=$#line;$j++){
	    print OUT   "  $line[$j]";
	    print       "  $line[$j]";
	    $ntatoms+=$line[$j];
	    $nsum[$j]=$ntatoms;
 	}
	print OUT   "\n";
	print       "\n";
	next LINE;
    }
    if($line[0] =~ /^[sS]/) {
        print OUT   $line."\n";
        print       $line."\n";
        $tag_sel = "True";
        next;}
    if($line[0] =~ /^[dD]/) {
        print OUT   "Cartesian\n";
        print       "Cartesian\n";
        next;}
    if($tag_sel == 'True'){
        $ncoord=6;
    }else{
        $ncoord=3;
    }
    for($j=0;$j<$ncoord;$j++){ $coord[$iatom][$j]=$line[$j];
        if($tag_sel == "True"){
            
        }
    }
    $iatom++;
    if($iatom >= $ntatoms){ last;}		
} continue{
    $i++;
}
close(IN); 

print "total number of atom :: $iatom\n";


for($i=0;$i<$ntatoms;$i++){
    ($x,$y,$z)=&xyz($coord[$i][0],$coord[$i][1],$coord[$i][2]);
    printf      "  %10.5f  %10.5f  %10.5f",  $x*$scale, $y*$scale, $z*$scale;
    printf OUT  "  %10.5f  %10.5f  %10.5f",  $x*$scale, $y*$scale, $z*$scale;
    if($tag_sel == "True"){
        printf      "  %2s  %2s  %2s", $coord[$i][3], $coord[$i][4], $coord[$i][5];
        printf OUT  "  %2s  %2s  %2s", $coord[$i][3], $coord[$i][4], $coord[$i][5];
    }
    print       "\n";
    print OUT   "\n";
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
