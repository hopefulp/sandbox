#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

##### default
%rem=(
    JOBTYPE => 'sp',
    EXCHANGE => 'hf',
    CORRELATION => 'rimp2',
    BASIS => 'cc-pvtz',
    SCF_MAX_CYCLES => '200',
    MOLDEN_FORMAT => 'TRUE'
);


##### IF THERE IS ARGUMENT, READ VARIABLE WHICH HAS "="
if($#ARGV >= 0){
    #print $#ARGV;
    for($i=0;$i<=$#ARGV;$i++){
	@iarg=split(/=/,$ARGV[$i]);
	if($#iarg != 1){
	    next;
	}
	if($iarg[0]="jobname"){
	    $jobname=$iarg[1];
	}else{	$rem{$iarg[0]}=$iarg[1];
	}
    }
}

##### if correlation is rimp2, add below
if($rem{'CORRELATION'}='rimp2'){ 
    $rem{'PURECART'}='1111';
    if($rem{'BASIS'}='cc-pvtz'){
    	$rem{'AUX_BASIS'}='rimp2-cc-pvtz';
    }
}

$newjob=$jobname;
print_hash();

# DFT - XYGJOS
#%rem=(
#    JOBTYPE => 'sp',
#    EXCHANGE => 'xygjos',
#    CORRELATION => '',
#    BASIS => 'aug-cc-pvtz',
#    AUX_BASIS => 'rimp2-cc-pvtz',
#    SCF_MAX_CYCLES => '200',
#    MOLDEN_FORMAT => 'TRUE'
#);

#operate on all the files in "dir.txt"
open(INd,"<dir.txt");
while($filename=<INd>){
    chomp($filename);
    $fname=$filename;
#$fname=$ARGV[0];
@fname=split(/\./,$fname);
$fname_pre=$fname[0];

$fin=$fname_pre."_$newjob.inp";
$fout=$fname_pre."_$newjob.out";
open(IN,"<$fname");
open(OUT,">$fin");

$print_tag="OFF";
@line_list=();
$i=0;
### save from $molecule to the EOF 
while($line=<IN>){
    @field=split(/\s+/,$line);
    if($field[0] eq "") {shift(@field);}
    if($field[0] eq "\$molecule"){
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

}

close(INd);

sub print_hash{
    my($key, $value);

    while(($key,$value)=each(%rem)){
	print "$key => $value\n";
    }
}
