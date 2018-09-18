#!/usr/bin/perl

$fin=$ARGV[0];

%atom_pairs=( 1 => '32', 2 => '29', 3 => '33', 4 => '25', 5 => '27', 6 => '36');

open(IN,"<$fin");

@first_coord=();
@second_coord=();
$i=0;
$iatom=0;
while($line=<IN>){

    chomp($line);
    @line=split(/\s+/, $line);
    if($line[0] eq ""){ shift(@line);}
    if($i==0 and $line[0] =~ /\D/){
	print "Error: In xyz format, first line is expected to be a number!\n";
	next;
    }
    next if($i<2);
    
    $iatom++;
    
    for $first ( keys %atom_pairs ){
	if($iatom == $first){
	    push(@first_coord,$line[1]);push(@first_coord,$line[2]);push(@first_coord,$line[3]);
	}
    }
    for $second ( values %atom_pairs ){
	if($iatom == $second){
	    # for second atom find first atom in the pair
	    for $first ( keys %atom_pairs ){
		if($atom_pairs{$first} == $iatom){
		    $pivot=$first;
		    last;
		}
	    }
	    for($k=0;$k<3;$k++){
		$second_coord[3*($pivot-1)+$k]=$line[$k+1];
	    }
    	}
    }
#    print join(" ",@second_coord)."\n";

} continue {
    $i++;
}

close(IN);

@distance=();
$dist_sum=0;
for($i=0;$i<=$#first_coord;$i+=3){
    $xdist=$first_coord[$i]-$second_coord[$i];
    $ydist=$first_coord[$i+1]-$second_coord[$i+1];
    $zdist=$first_coord[$i+2]-$second_coord[$i+2];

    $sqdist=$xdist*$xdist+$ydist*$ydist+$zdist*$zdist;
    $dist=sqrt($sqdist);
    $dist_sum+=$dist;
    push(@distance,$dist);
    $turn=int($i/3);
#    print "$turn -th distance is $dist\n";
    
}

$ave_dist=$dist_sum/($#distance+1);
print "Ave distance is $ave_dist\n";


