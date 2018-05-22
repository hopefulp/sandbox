#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.inp

$atom=$ARGV[1];

##### default
%rem=(
    JOBTYPE => 'sp',
    EXCHANGE => 'hf',
    CORRELATION => 'rimp2',
    BASIS => 'gen',
    SCF_MAX_CYCLES => '200',
    MOLDEN_FORMAT => 'TRUE',
    PURECART => '1111',
    AUX_BASIS => 'rimp2-cc-pvtz'
);

print_hash();

system("get_basis.pl $atom > $atom.bas");

@basis=();
open(IN,"<$atom.bas");
while($line=<IN>){
     push(@basis,$line);
}

$filename=$ARGV[0];
open(OUT,">>$filename");

print OUT "\n\$basis\n";
print OUT @basis;
print OUT "\$end\n";
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
