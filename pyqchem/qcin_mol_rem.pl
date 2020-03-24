#!/usr/bin/perl
# written by Joonho Park
# read a.out and write a.mol, read the initial geometry
if($#ARGV < 0){
    print "scratch qchem input file from a.out \n";
    print "Usage:: $0 r=remfile m=geofile i=outfile\n";
    print "read input information\n";
    exit(0);
}

for($i=0;$i<=$#ARGV;$i++){
    @iarg=split(/=/,$ARGV[$i]);
    if($iarg[0] =~ /r/){
    	$remfile=$iarg[1];
    }elsif($iarg[0] =~ /i/){
	    $inpfile=$iarg[1];
    }elsif($iarg[0] =~ /m/){
    	$geofile=$iarg[1];
    }else{
        $geofile=$ARGV[$i];
        #$inpfile=$ARGV[$i];
    }
}

if ($remfile eq ""  and $geofile ne ""){
    @filename=split(/\./, $geofile);
    $inpfile=$filename[0].".mol";
}elsif ($remfile ne "" and $geofile eq ""){
    @filename=split(/\./, $remfile);
    $inpfile="rem.".$filename[0];
}elsif ($remfile ne "" and $geofile ne ""){
    @filename=split(/\./, $geofile);
    $inpfile=$filename[0].".in";
}elsif($inpfile ne ""){
    @filename=split(/\./, $inpfile);
    $remfile=$inpfile;
    $geofile=$inpfile;
    $inpfile=$filename[0].".in";
}else{
    print "Error in input\n";
    exit(1);
}    

#print "geo file = $geofile, outfile = $inpfile\n";
open(OUT,">$inpfile");
print "outfile is $inpfile\n";
$flag="OFF";
@bas_line=();
@input=($remfile, $geofile);
@keyword=("rem", "molecule");
#print @keyword,"\n";
for($i=0;$i<2;$i++){
    open(IN,"<$input[$i]");
    while($line=<IN>){
    	#@field=split(/\s+/,$line);
	    if($flag eq "OFF"){			# wait until signal turns on
    	    if($line =~ /$keyword[$i]/){     #/$keyword[$i]/){
	    	    $flag="ON";
	    	    push(@bas_line,$line);
	        }
            next;
        }else{
            if($line =~/\$end/){
                push(@bas_line,$line);
                $flag="OFF"; 
		        last;
            }else{
                #### remove blank line if it is in the $molecule block
                if($line =~ /^$/) { next; }
                push(@bas_line,$line);
            }
        }
    }
    close(IN);
    if($i==$#keyword+1){last;}
}
print @bas_line;
print OUT @bas_line;

close(OUT); 
