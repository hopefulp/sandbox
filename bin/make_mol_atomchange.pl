#!/usr/bin/perl
# written by Joonho Park
# read  a.dat and a.inp then write me.mol
# a.inp includes mol coordinates

$fdat=$ARGV[0];
$fin=$ARGV[1];
$fmodel=$ARGV[2];
$charge=-1;

if($#ARGV < 2){
    print "argument 1st dat list, 2nd mol file, 3rd model name\n";
    exit(0);
}

open(IN,"<$fdat");
$find="NO";
$flag="OFF";
$key="\$molecule";
$key2="\$end";
### read three atoms's coordinates
while($line=<IN>){
#    $i++;
    #print $line;
    @field=split(/\s+/,$line);
    if($field[0] eq "") {shift(@field);}
    $atom=$field[0];
    $multip=$field[1];
    print "Atom:$atom\tMultip:$multip\n";
    $fout="${atom}_$fmodel.mol";
    open(OUT,">$fout");
    open(IN2,"<$fin");

    @f_line=();
    $i=0;
    while($row=<IN2>){
	@part=split(/\s+/,$row);
	if($part[0] eq "") {shift(@part);}
        if($flag eq "OFF"){			# wait until signal turns on
            if($part[0] eq $key){
##   	        print $i."\n";
                $flag="ON";
	        push(@f_line,$row);
            }else{ next;}
        }else{
            $i++;
            if($part[0] eq $key2){
                $flag="OFF"; 
		push(@f_line,$row);
		next;
            }elsif($i<2){
		$new_line="\t$charge\t$multip\n";
		push(@f_line,$new_line);
                next;
            }elsif($i==2){
		$new_line="  $atom\t$part[1]\t$part[2]\t$part[3]\n";
                push(@f_line,$new_line);
	    }else{ push(@f_line,$row);}
        }
    }
    close(IN2);
    print OUT @f_line;
    close(OUT);
}
close(IN); 
