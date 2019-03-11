package IOFile2;

use Exporter;

@ISA=qw(Exporter);
@EXPORT=qw(file1,file2);

sub file1{
    my($fin)=@_;
    my($fout,@name);

    if($fin eq ""){
	print "Error: input file\n";
	exit;
    }
#    if($suffix eq ""){
#	$suffix="suff";
#	print "Suffix is added by $suffix\n";
#    }

    if($fin =~ /\./){
    	@name=split(/\./,$fin);
    	$fout=$name[0].".$suffix";
    }else{
    	$fout=$fin.".$suffix";
    }
   
    print "IOFILE::OUTFILE is $fout\n";

    return $fout;
}

sub file2{
    my($fin)=@_;
    my($fout,@name);

    if($fin eq ""){
	print "Error: input file\n";
	exit;
    }

    @name=split(/\./,$fin);
    if($name[1] eq ""){
        $name[1]="mol";
    }
    
    if($name[1] eq "mol"){
        $suffix = "xyz";
    }elsif($name[1] eq "xyz"){
        $suffix = "mol";
    }else{
        print "I cannot understand input file format: function mol2xyz or xyz2mol\n";
    }

    $fout=$name[0].".$suffix";
   
    print "IOFILE::OUTFILE is $fout\n";

    return $fout, $suffix;
}

sub file3{
    my($fin)=@_;
    my($fout,@name);

    if($fin eq ""){
        print "Error: input file\n";
        exit;
    }

    @name=split(/\./,$fin);

    if($name[1] eq "Amol"){
        $suffix = "mol";
    }else{
        print "input file is supposed to be Accelrys mol file\n";
	exit (1);
    }

    #$fout=$name[0].".$suffix";
    $fout=$name[0];

    print "IOFILE::OUTFILE is $fout\n";

    return $fout, $suffix;
}




1;
