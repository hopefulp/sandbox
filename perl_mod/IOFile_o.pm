package IOFile;

use Exporter;

@ISA=qw(Exporter);
@EXPORT=qw(file1, Usage);

sub file1{
    my($fin,$suffix)=@_;
    my($fout,@name);

    if($fin eq ""){
	print "Error: input file\n";
	exit;
    }
    if($suffix eq ""){
	$suffix="suff";
	print "Suffix is added by $suffix\n";
    }

    if($fin =~ /\./){
    	@name=split(/\./,$fin);
    	$fout=$name[0].".$suffix";
    }else{
    	$fout=$fin.".$suffix";
    }
   
    print "IOFILE::OUTFILE is $fout\n";

    return $fout;
}

sub Usage{
    my($app, @args)=@_;
    my($i, $nargs);

    print "Usage::$app ";
    for($i=0; $i<=$#args; $i++){
        print "$args[$i] ";
    }
    print "\n";
    exit(1);
}

1;
