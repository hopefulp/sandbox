#!/usr/bin/perl
# written by Joonho Park
# read INCAR and overwrite INCAR

print "Warning: Overwrite INCAR\n";

$job=$ARGV[0];
$fout="tmp.incar";

if($job =~ /=/){
    @jobs=split(/=/, $job);
    $job=$jobs[0];
    $job_method=$jobs[1];
}

open(IN,"< INCAR");
open(OUT,">$fout");

@f_line=();
### read three atoms's coordinates
$i=0;

$key1="GGA";
$key2="LDAUU";
LINE: while($line=<IN>){
#    $i++;
    #print $line;
    @line=split(/\s+/,$line);
    if($line[0] eq "") { shift @line; }

    if($line[0] =~ /$key1/i){
	if($job eq $l_jobs[1]){
	    print OUT "$key1 = $job_method\n";
	    next LINE;
	}
    }
    if($line =~ /$key2/i){
	if($job eq "U"){
	    print OUT "$key2 = $job_method 0 0 0\n";
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

