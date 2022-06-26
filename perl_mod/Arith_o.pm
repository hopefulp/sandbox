# by J. Park

package Arith;

use Exporter;
use List::Util qw( min max );

@ISA=qw(Exporter);
@EXPORT=qw(dist, min_list, dist_lattice, dist_lattice_z, dist_cubic, Gaussian_plot);


sub dist{
    my($x1,$y1,$z1,$x2,$y2,$z2)=@_;
    my($dist_sq, $dist);

    if($#_ != 5){
	print "Error: The number of arguments is not 6 for 3D distance\n";
	exit;
    }

    $dist_sq=($x2-$x1)*($x2-$x1)+($y2-$y1)*($y2-$y1)+($z2-$z1)*($z2-$z1);
    $dist=sqrt($dist_sq);

#    print "distance = ", $dist, "\n";

    return($dist);
}

### This works for 1d on z-direction
sub dist_lattice_z {
    my($lattice_a,$lattice_b,$lattice_c,$atom1,$atom2)=@_;
    my(@a_latt,@b_latt,@c_latt);
    my(@acos,@bcos,@ccos,$a_length,$b_length,$c_length,$a_half_leng,$b_half_leng,$c_half_leng);
    my(@a1_coord,@a2_coord,@replica_coord);
    my($d_x,$d_y,$d_z,$d_x2,$d_y2,$d_z2,$dist,$dist_crit);
    my($dist0,$dist_x1,$dist_y1,$dist_z1,$tag);

    $dist_crit=3;

    @a_latt=@{$lattice_a};    @b_latt=@{$lattice_b};	@c_latt=@{$lattice_c};
    $a_length=sqrt($a_latt[0]*$a_latt[0]+$a_latt[1]*$a_latt[1]+$a_latt[2]*$a_latt[2]);
    $b_length=sqrt($b_latt[0]*$b_latt[0]+$b_latt[1]*$b_latt[1]+$b_latt[2]*$b_latt[2]);
    $c_length=sqrt($c_latt[0]*$c_latt[0]+$c_latt[1]*$c_latt[1]+$c_latt[2]*$c_latt[2]);
    $a_half_leng=$a_length/2; $b_half_leng=$b_length/2; $c_half_leng=$c_length/2;
#    $acos[0]=$a_latt[0]/$a_length;  $bcos[0]=$b_latt[0]/$b_length;   $ccos[0]=$c_latt[0]/$c_length;
#    $acos[1]=$a_latt[1]/$a_length;  $bcos[1]=$b_latt[1]/$b_length;   $ccos[1]=$c_latt[1]/$c_length;
#    $acos[2]=$a_latt[2]/$a_length;  $bcos[2]=$b_latt[2]/$b_length;   $ccos[2]=$c_latt[2]/$c_length;

    @a1_coord=@{$atom1}; @a2_coord=@{$atom2};

    #print "sub:: ",join(" ",@a1_coord)," \n";
    #print "sub:: ",join(" ",@a2_coord)," \n";
    #print "sub1: $a_length $b_length $c_length\n";
    ###### IF the orgin is (0,0,0), you can just add lattice vector.
    ### +Z
    $d_x=&dist_1d($a1_coord[0],$a2_coord[0]);
    $d_y=&dist_1d($a1_coord[1],$a2_coord[1]);
    $d_z=&dist_1d($a1_coord[2],$a2_coord[2]);
    $dist=sqrt($d_x*$d_x+$d_y*$d_y+$d_z*$d_z);
    #print "sub dist: $d_x $d_y $d_z $dist\n"; 
    $tag=0;
#    if($dist >= 4){
    ### d_z should be projected to c-direction but Now c is similar to z
    ### atom1 is reference and move atom2 to the close to atom1
    	#print "sub2: $d_z  $c_half_leng\n";
    	if($d_z >= $c_half_leng){
	    if($a1_coord[2] < $c_half_leng){	# move atom2 to -z
	    	$replica_coord[0]=$a2_coord[0]-$c_latt[0];
	    	$replica_coord[1]=$a2_coord[1]-$c_latt[1];
	    	$replica_coord[2]=$a2_coord[2]-$c_latt[2];
		$tag=-1;
	    }else{				# move atom2 to +z
	    	$replica_coord[0]=$a2_coord[0]+$c_latt[0];
	        $replica_coord[1]=$a2_coord[1]+$c_latt[1];
	        $replica_coord[2]=$a2_coord[2]+$c_latt[2];
		$tag=1;
	    }
	    $dist=&dist_3d(\@a1_coord,\@replica_coord);
	    if($dist > $dist_crit){ print "Error with long distance::&dist_lattice_z\n"; exit(112); }
	    ### the value in main program should be changed but why?
	    @{$atom2}=@replica_coord;
	}
   return($tag, $dist);
}

sub dist_lattice {
    my($lattice_a,$lattice_b,$lattice_c,$atom1,$atom2)=@_;
    my(@a_latt,@b_latt,@c_latt);
    my(@acos,@bcos,@ccos,$a_length,$b_length,$c_length,$a_half_leng,$b_half_leng,$c_half_leng);
    my(@a1_coord,@a2_coord,@replica_coord);
    my($d_x,$d_y,$d_z,$d_x2,$d_y2,$d_z2,$dist,$dist_crit,$dist1,$dist2);
    my($dist0,$dist_x1,$dist_y1,$dist_z1,$tag);

    $dist_crit=3;

    @a_latt=@{$lattice_a};    @b_latt=@{$lattice_b};	@c_latt=@{$lattice_c};
    $a_length=sqrt($a_latt[0]*$a_latt[0]+$a_latt[1]*$a_latt[1]+$a_latt[2]*$a_latt[2]);
    $b_length=sqrt($b_latt[0]*$b_latt[0]+$b_latt[1]*$b_latt[1]+$b_latt[2]*$b_latt[2]);
    $c_length=sqrt($c_latt[0]*$c_latt[0]+$c_latt[1]*$c_latt[1]+$c_latt[2]*$c_latt[2]);
    $a_half_leng=$a_length/2; $b_half_leng=$b_length/2; $c_half_leng=$c_length/2;
#    $acos[0]=$a_latt[0]/$a_length;  $bcos[0]=$b_latt[0]/$b_length;   $ccos[0]=$c_latt[0]/$c_length;
#    $acos[1]=$a_latt[1]/$a_length;  $bcos[1]=$b_latt[1]/$b_length;   $ccos[1]=$c_latt[1]/$c_length;
#    $acos[2]=$a_latt[2]/$a_length;  $bcos[2]=$b_latt[2]/$b_length;   $ccos[2]=$c_latt[2]/$c_length;

    @a1_coord=@{$atom1}; @a2_coord=@{$atom2};

    #print "sub:: ",join(" ",@a1_coord)," \n";
    #print "sub:: ",join(" ",@a2_coord)," \n";
    #print "sub1: $a_length $b_length $c_length\n";
    ###### IF the orgin is (0,0,0), you can just add lattice vector.
    ### +Z
    $d_x=&dist_1d($a1_coord[0],$a2_coord[0]);
    $d_y=&dist_1d($a1_coord[1],$a2_coord[1]);
    $d_z=&dist_1d($a1_coord[2],$a2_coord[2]);
    $dist1=sqrt($d_x*$d_x+$d_y*$d_y+$d_z*$d_z);
    #print "sub dist1: $d_x $d_y $d_z $dist1\n"; 
    $tag=0;
#    if($dist1 >= 4){
    ### d_z should be projected to c-direction but Now c is similar to z
    ### atom1 is reference and move atom2 to the close to atom1
    	#print "sub2: $d_z  $c_half_leng\n";
    	if($d_z >= $c_half_leng){
	    if($a1_coord[2] < $c_half_leng){	# move atom2 to -z
	    	$replica_coord[0]=$a2_coord[0]-$c_latt[0];
	    	$replica_coord[1]=$a2_coord[1]-$c_latt[1];
	    	$replica_coord[2]=$a2_coord[2]-$c_latt[2];
		$tag=-1;
	    }else{				# move atom2 to +z
	    	$replica_coord[0]=$a2_coord[0]+$c_latt[0];
	        $replica_coord[1]=$a2_coord[1]+$c_latt[1];
	        $replica_coord[2]=$a2_coord[2]+$c_latt[2];
		$tag=1;
	    }
	    $dist2=&dist_3d(\@a1_coord,\@replica_coord);
	    if($dist2>$dist1){  $dist=$dist1;}
	    else{ 		$dist=$dist2;}
	    ### the value in main program should be changed
	    if($tag != 0 && $dist < $dist_crit){
		print "2nd atom coordinate was change\n";
		print Dumper $atom2;
		print Dumper \@replica_coord;
	    	@{$atom2}=@replica_coord;
	    }
	}else{$dist=$dist1;}
   return($tag, $dist);
}

sub dist_cubic {
    # arguments are references 
    my($lattice,$atom1,$atom2)=@_;
    my($d_x,$d_y,$d_z,$dist);
    my(@a1_coord, @a2_coord, @lattice);

    @a1_coord=@{$atom1}; @a2_coord=@{$atom2}; @lattice=@{$lattice};


    #print "sub1: $a_length $b_length $c_length\n";
    ###### IF the orgin is (0,0,0), you can just add lattice vector.
    ### +Z
    $d_x=&dist_1d_lat($a1_coord[0],$a2_coord[0],$lattice[0]);
    $d_y=&dist_1d_lat($a1_coord[1],$a2_coord[1],$lattice[1]);
    $d_z=&dist_1d_lat($a1_coord[2],$a2_coord[2],$lattice[2]);
    $dist=sqrt($d_x*$d_x+$d_y*$d_y+$d_z*$d_z);
    #print "sub dist1: $d_x $d_y $d_z $dist1\n"; 
   return($dist);
}

sub dist_1d_lat {
    my ($a, $b, $lattice)=@_;
    my ($dist,$a_image, @dist, $h_latt);

    $h_latt= $lattice/2.;
    $dist=abs($a-$b);
    if($dist < $h_latt){
	    return $dist;
    }else{
	    push (@dist, $dist);
        $a_image=$a+$lattice;
        $dist=abs($a_image-$b);
        push(@dist, $dist);
        $a_image=$a-$lattice;
        $dist=abs($a_image-$b);
        push(@dist, $dist);
    }
    #print join(" ", @dist),"\n";
    $min = min @dist;
    return $min;
}

sub dist_3d {
    my ($a, $b)=@_;
    my (@a, @b, $dist);

    @a=@{$a};
    @b=@{$b};

    $dist=sqrt(($a[0]-$b[0])*($a[0]-$b[0])+($a[1]-$b[1])*($a[1]-$b[1])+($a[2]-$b[2])*($a[2]-$b[2]));
    return $dist;
}
 
sub dist_1d {
    my ($a, $b)=@_;
    my ($dist);

    $dist=abs($a-$b);
    return $dist;
}

sub min_list {
    my(@list)=@_;
    my($i,$min);

    $min=1e10;
    for($i=0;$i<=$#list;$i++){
	if($list[$i]<$min){
	    $min=$list[$i];
    	}
    }
    return $min;
}

1;
