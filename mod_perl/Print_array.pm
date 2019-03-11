package Print_array;

use Exporter;

@ISA=qw(Exporter);
@EXPORT=qw(compare_two cal_dist print_hash print_list print_matrix);

sub compare_two {
    my($i,$j,@coord)=@_;
    my($k);

    for($k=0;$k<3;$k++){
        print " comp: $coord[$i][$k] $coord[$j][$k]\n"
    } print "\n";
}

sub  cal_dist {
    my($i,$j,@coord)=@_;
    my($k,$sum2,$dist);

    $sum2=0;
    for($k=0;$k<3;$k++){
        $sum2+=($coord[$i][$k]-$coord[$j][$k])*($coord[$i][$k]-$coord[$j][$k]);
    }
    $dist=sqrt($sum2);
}

#sub print_hash {
#    my($key, $value);
#
#    while(($key,$value) = each @_;
#    	print "key $key, value $value\n";
#}
sub print_list {
    my(@list)=@_;
    my($i);

    print join(" ",@list)."\n";
}

sub print_matrix {
    my($nrow,$ncol, @matrix)=@_;
    my($i,$j);

    for($i=0;$i<$nrow;$i++){
        for($j=0;$j<$ncol;$j++){
            print "  $matrix[$i][$j]";
        }
        print "\n";
    }
}


1;
