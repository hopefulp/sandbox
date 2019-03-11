### written by J. Park

package GetHash;
use Exporter;

@ISA=qw(Exporter);
@EXPORT=qw(getvalue get_2values);

### from the input such as "l=3", split l and 3
sub getvalue{
    my($arg,$delimiter)=@_;
    my($key,$value);

    if($arg eq ""){
        print "Use module UseHash::getvalue: &UseHash::getvalue(\"ab=number\",\"$delimiter\") \n";
        exit;
    }

    ($key,$value)=split(/$delimiter/,$arg);
#    ($key,$value)=split(/=/,$arg);
#    print $key,"del", $value,"\n";   
 
    return ($key, $value);

}


### from the input such as "a:3:4", split a, 3, and 4 then return (3, 4)
sub get_2values{
    my($arg,$delimiter)=@_;
    my(@list,$prefix,$start,$end);

    if($arg eq ""){
        print "Use module UseHash::get_2values: &UseHash::get_2values(\"a:n_start:n_end\",\"$delimiter\") \n";
        exit;
    }

    @list=split(/:/,$arg);

    if($list[0] !~ /\d/){
	$prefix=shift(@list);
    }
    $start=$list[0];
    $end=$list[1];
    if($end eq ""){
	$end = $start;
    }
    print "prefix: ",$prefix,"  start: ",$start,"  end: ", $end,"\n";

    return ($start,$end);
}

jp;
