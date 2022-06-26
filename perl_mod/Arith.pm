# by J. Park

package Arith;

use Exporter;
use List::Util qw( min max );
use Data::Dumper qw(Dumper);


@ISA=qw(Exporter);
@EXPORT=qw(dist, dist_cubic, Gaussian_plot);

$PI=3.14159;

### receive 1-D data and plot
sub Gaussian_plot {
    my($sigma, $min, $max, $np, @ndata)=@_;
    my(@f,@f2, @tmp, $sum, $ave, $intensity, $prob );
    my($i, $x_i, $x_f, $dx, $xi);

    $sum=0;
#    $np=200;
#    $min = min @ndata;
#    $max = max @ndata;
#    for($i=0; $i<=$#points; $i++){
#        $sum+=$points[$i];
#    }
#    $ave=$sum/$#points;
    $x_i=$min; #-5*$sigma;
    $x_f=$max; #+5*$sigma;
    $dx=($x_f-$x_i)/$np;

    @f=();
    @f2=();
    @tmp=();
    for($i=0;$i<$np;$i++){
        $xi=$x_i+$i*$dx;
        $intensity=0.;
        for($j=0;$j<=$#ndata;$j++){
            $intensity+=f_gaussian($xi,$ndata[$j], $sigma);
        }
        $prob=$intensity/$#ndata;
        push(@f,$prob);
        #print "$xi $intensity\n";
        @tmp=($xi, $prob, $intensity);
        #print @tmp, "\n";
        push(@{$f2[$i]}, @tmp);
    }        
    #print join("\n",@f), "\n";
    #print Dumper \@f2;

    return \@f2;
}


sub f_gaussian{
    my ($x,$mu, $sigma)=@_;
    my $f;

    if( $mu-3*$sigma < $x and $x < $mu+3*$sigma){
        $f=1./($sigma*sqrt(2.*$PI))*exp(-0.5*($x-$mu)*($x-$mu)/($sigma*$sigma));
        return $f;
    }else{  return 0;}
}

   


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

1;
