#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

$filename=$ARGV[0];
@fname=split(/\./,$filename);
$fname_pre=$fname[0];

$newjob="bas";
$fin=$fname_pre."_$newjob.inp";		#qchem input file
$fout=$fname_pre."_$newjob.out";	#qchem output file
open(IN,"<$filename");
open(OUT,">$fin");

$print_tag="OFF";
@line_list=();
$i=0;
### save from $molecule to the EOF 
while($line=<IN>){
    @field=split(/\s+/,$line);
    if($field[0] eq "") {shift(@field);}
    if($field[0] eq "\$basis"){
	$print_tag="ON";
    }

    if($print_tag eq "ON"){
	push(@line_list, $line);
	#print $line."print $i\n";
    }
    $i++;
}
close(IN);
print OUT  "\$rem\n";
while(($key,$value) = each(%rem)){
    print OUT "    $key	$value\n";
}
print OUT "\$end\n";
for($i=0;$i<$#line_list;$i++){
    print OUT $line_list[$i];
}
close(OUT);
#system("qchem $fin $fout &"); 

################## sub routine
sub print_hash{
    my($key, $value);

    while(($key,$value)=each(%rem)){
	print "$key => $value\n";
    }
}
