#!/usr/bin/perl
# written by Joonho Park
# read INCAR and overwrite INCAR

@l_jobs=qw(cont GGA); #job list

if($#ARGV < 0){
    print "Usage:: $0  job_type\n";
    print " jobs:: ",join(" ", @l_jobs),"\n";
    print "     cont:: if cont => ISTART 1 ICHARG 0 LAECHG .TRUE.\n";
    print "	       not  =>        0        2        .FALSE.\n";
    print "     gga :: input method; PE RP etc\n";	
    exit(1);
}

#use List::MoreUtils 'any';

print "Warning: Overwrite INCAR\n";

$job=$ARGV[0];
$mag=$ARGV[1];
$fout="tmp.incar";

if($job =~ /=/){
    @jobs=split(/=/, $job);
    $job=$jobs[0];
    $job_method=$jobs[1];
}


for($i=0;$i<=$#l_jobs;$i++){
    if($job eq $l_jobs[$i]){
	last;
    }
    if($i==$#l_jobs){
	print "There is no job for $job\n";
	exit(2);
    }
}

open(IN,"< INCAR");
open(OUT,">$fout");

@f_line=();
### read three atoms's coordinates
$i=0;

$key1="ISTART";
$key2="ICHARG";
$key3="GGA";
$key4="LAECHG";
LINE: while($line=<IN>){
#    $i++;
    #print $line;
    @line=split(/\s+/,$line);
    if($line[0] eq "") { shift @line; }

    if($line[0] =~ /$key1/i){
	if($job eq "cont"){
		print OUT "$key1 = 1\n";
	}else{  print OUT "$key1 = 0\n";}
	next LINE;
    }
    if($line[0] =~ /$key2/i){
	if($job eq "cont"){
                print OUT "$key2 = 0\n";
        }else{  print OUT "$key2 = 2\n";}
	next LINE;
    }
    if($line[0] =~ /$key3/i){
	if($job eq $l_jobs[1]){
	    print OUT "$key3 = $job_method\n";
	    next LINE;
	}
    }
    if($line =~ /$key4/i){
	if($job eq "cont"){
	    print OUT "$key4 = .TRUE.\n";
	    next LINE;
	}
    }
    print OUT $line;
    next LINE;
}
continue{
    $i++;
}
close(IN); 
close(OUT);

system("cp tmp.incar INCAR");

