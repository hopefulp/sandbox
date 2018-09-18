#!/usr/bin/perl
# written by Joonho Park
# read OUTCAR and print magnetization

if($#ARGV < 2){
    print "Usage:: $0 dir 1st_atom last_atom\n";
    exit(1);
}

$fin="OUTCAR";
$fout="magnetization";

$dir=$ARGV[0];
$iatom=$ARGV[1];
$fatom=$ARGV[2];
$mag= $ARGV[3];
$half_atom=($fatom-$iatom+1)/2;

open(IN,"<$dir/$fin");
open(OUT,">$fout");

$key_word="magnetization";

$tag_p="NO";
$tag_mag_p="YES";
if($mag eq ""){
    $mag = "AFM";
    print "Mag: $mag for default\n";
}

if($#ARGV >= 3){
    $tag_mag_p="YES";
}


$n_atom=$fatom-$iatom+1;
$tag="NO";
$i=0;
$i_tag=0;
@mag=();
### read OUTCAR file
while($line=<IN>){
    chomp($line);
    if($tag eq "NO" and $line =~ m/magnetization/){
    	@field=split(/\s+/,$line);
    	if($field[0] eq "") {shift(@field);}
	if($field[0] eq $key_word){
	    if($tag_p eq "YES") {print "$line\n";}
	    if($tag_p eq "YES") {print "start read\n";}
	    $tag="YES";
	    next;
	}
    }elsif($tag eq "YES"){
	if($i_tag <= 2){
	    $i_tag++;
	    next;
	}else{
	    @field=split(/\s+/,$line);
	    if($field[0] eq "") {shift(@field);}
	    if($iatom <= $field[0] and $field[0] <= $fatom){
		#print "$i: $tag,  $field[0] $field[4]\n";
		if($tag_mag_p eq "YES") {printf "%2d\t%10.3f\n", $field[0], $field[4];}
		push(@mag, $field[4]);
	    }
	    if($field[0] == $fatom){
		$tag="NO";
		last;
	    }
	}
   }
}continue{
    $i++;
}

close(IN); 
close(OUT);



if($tag_mag_p eq "YES") {print join(" ",@mag),"\n";}
$sum_mag=0;
#print "MAG: $mag\n";
for($i=0;$i<$n_atom;$i++){
    $sum_mag+=abs($mag[$i]);
    if($i<$half_atom){
	if($mag eq "AFM"){
	    if($i%2==0 && $mag[$i]<0) {print "Error in $mag\n"; exit(99);}
	    elsif($i%2==1 && $mag[$i]>0) {print "Error in $mag\n"; exit(99);}
	}else{
	    if($mag[$i]<0) {print "Error in FM\n"; exit(99);}
	}
    }else{
	if($mag eq "AFM"){
	    if($i%2==0 && $mag[$i]>0) {print "Error in $mag\n"; exit(99);}
	    elsif($i%2==1 && $mag[$i]<0) {print "Error in $mag\n"; exit(99);}
	}else{
	    if($mag[$i]>0) {print "Error in FM\n"; exit(99);}
	}
    }
}

$ave_mag=$sum_mag/$n_atom;
printf "Ave magnetization %10.3f\n", $ave_mag;


