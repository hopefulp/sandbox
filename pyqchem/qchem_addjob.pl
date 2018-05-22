#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.inp

##### default
%rem=(
    JOBTYPE => 'sp',
    EXCHANGE => 'hf',
    CORRELATION => 'rimp2',
    BASIS => 'cc-pvtz',
    SCF_MAX_CYCLES => '200',
    MOLDEN_FORMAT => 'TRUE',
    PURECART => '1111',
    AUX_BASIS => 'rimp2-cc-pvtz'
);


print_hash();


$filename=$ARGV[0];
$fin=$filename.".tmp";			#qchem temparary input file
open(OUT,">>$filename");

print OUT "\n@@@\n\n";
print OUT "\$molecule\n    read\n\$end\n";
print OUT "\$rem\n";
while(($key,$value) = each(%rem)){
    print OUT "    $key	$value\n";
}
print OUT "\$end\n";

close(OUT);

################## sub routine
sub print_hash{
    my($key, $value);

    while(($key,$value)=each(%rem)){
	print "$key => $value\n";
    }
}
