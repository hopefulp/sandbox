#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a_newjob.input and run qchem a_newjob.inp a_newjob.out
# if $molecule ne "TRUE", read "mol.inp" in the working directory
# if basis eq gen and the basis is not on the list of if-sentence make atom.bas file


################################# DEFAULT REM
%rem=(
    JOBTYPE => 'sp',
    EXCHANGE => 'hf',
    CORRELATION => 'rimp2',
    BASIS => 'cc-pvtz',
    SCF_MAX_CYCLES => '200',
    MOLDEN_FORMAT => 'TRUE',
);

$molecule="TRUE";

#################################### GET REM FROM ARGUMENTS
##### IF THERE IS ARGUMENT, READ VARIABLE WHICH HAS "="
##### read jobname if not, jobname = CORRELATION

if($#ARGV >= 0){
    #print $#ARGV;
    for($i=0;$i<=$#ARGV;$i++){
	@iarg=split(/=/,$ARGV[$i]);
	if($#iarg != 1){
	    $filename=$ARGV[$i];
	    next;
	}
	if($iarg[0] eq "jobname"){
	    $jobname=$iarg[1];
	}elsif($iarg[0] eq "molecule"){
	    $molecule=$iarg[1];
	}else{	
	    delete $rem{$iarg[0]} if exists $rem{$iarg[0]};
	    $rem{$iarg[0]}=$iarg[1];
	}
    }
}

if($rem{'JOBTYPE'} eq 'opt'){
    $rem{'GEOM_OPT_MAX_CYCLES'} = '200';
}

##### if correlation is rimp2, add below
if($rem{'CORRELATION'} eq 'rimp2'){ 
    $rem{'PURECART'}='1111';
    if($rem{'BASIS'} eq 'cc-pvtz'){
    	$rem{'AUX_BASIS'}='rimp2-cc-pvtz';
    }elsif($rem{'BASIS'} eq 'vdz(d)'){
	$rem{'AUX_BASIS'}='rimp2-VDZ';
    }
}
if($jobname eq ""){
    $newjob=$rem{'CORRELATION'}
}else {$newjob=$jobname;
}

#operate on all the files in "dir.txt"
#$filename=$ARGV[0];
@fname=split(/\./,$filename);
$fname_pre=$fname[0];

$fin=$fname_pre."_$newjob.inp";		#qchem input file
$fout=$fname_pre."_$newjob.out";	#qchem output file
open(OUT,">$fin");

$print_tag="OFF";
@line_list=();
$i=0;


# writing at the end

#################################### MOLECULE PART
########################## read input file and save molecule and atom species

$flag="OFF";
@atom_line=();
@atoms=();
    open(IN,"<$filename");
    while($line=<IN>){
        @field=split(/\s+/,$line);
	if($field[0] eq ""){shift(@field);}
        if($flag eq "OFF"){                     # wait until signal turns on
            if($field[0] eq "\$molecule"){
                $flag="ON";
                push(@atom_line,$line);
            }else{ next;}
        }else{
            if($field[0] eq "\$end"){
                push(@atom_line,$line);
                $flag="OFF";
                last;
            }else{
                push(@atom_line,$line);
		if($rem{'BASIS'} eq 'gen' or $rem{'BASIS'} eq 'GEN'){
		    $newatom="YES";
		    for($i=0;$i<=$#atoms;$i++){
			if($field[0] eq $atoms[$i]){
			    $newatom="NO";
			    last;
			}
		    }
		    if($newatom eq "YES"){ # if it is new and it is not number
			if($field[0]=~/\D/) {push(@atoms,$field[0]);}
		    }
		}
	    }
        }
    }
    close(IN);
print @atoms;
print "\n";

##### molecule from input file or read from outside when qchem runs
if($molecule eq "TRUE"){
    print OUT @atom_line;
}else{ 
    print OUT "\$molecule\n";
    print OUT "    read mol.inp\n";
    print OUT "\$end\n";
}

############################ ADD BASIS
# if atoms are not saved, this is skipped by for sentence

if($rem{'BASIS'} eq 'gen' or $rem{'BASIS'} eq 'GEN'){

$common_basis="cc-pvtz";
print OUT "\$basis\n";
for($i=0;$i<=$#atoms;$i++){
    $atom_list="ISNOT";
    ################### add basis in the qchem
    if($atoms[$i] eq "H" or $atoms[$i] eq "C" or $atoms[$i] eq "N" or $atoms[$i] eq "O" or \
	$atoms[$i] eq "S" ){
	print OUT "   $atoms[$i]    0\n";
	print OUT "   $common_basis\n";
	print OUT "   ****\n";
	$atom_list="IS"; next;
    }
    if($atom_list eq "ISNOT"){
	system("/home/joonho/bin/get_basis.pl $atoms[$i] 2 >/dev/null");
    }
    open(IN,"<$atoms[$i].bas");
    while($line=<IN>){
	print OUT $line;
    } 
    close(IN);
    if($rem{'CORRELATION'} eq "MP2" or $rem{'CORRELATION'} eq "mp2"){
	$rem{'MEM_STATIC'} = "1000";	# this default is effective for Mg
	if($atoms[$i] ne "Mg"){
	    delete $rem{"MEM_STATIC"} if exists $rem{"MEM_STATIC"};
	    $rem{'MEM_STATIC'} = "2000";
	    $rem{'MEM_TOTAL'}  = "3000";
	    if($atoms[$i] eq "Cu"){
	    	$rem{'MAX_CIS_CYCLES'} = "100";
	    }  
	}
    }
}
print OUT "\$end\n";

}

print OUT  "\$rem\n";
while(($key,$value) = each(%rem)){
    print OUT "    $key	$value\n";
}
print OUT "\$end\n";

close(OUT);

print_hash();

system("qchem $fin $fout &"); 

################## sub routine
sub print_hash{
    my($key, $value);

    while(($key,$value)=each(%rem)){
	print "$key => $value\n";
    }
}
