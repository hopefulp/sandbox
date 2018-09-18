#!/usr/bin/perl
# written by Joonho Park

$fname_pre=$ARGV[0];
$Aatom=$ARGV[1];
$Batom=$ARGV[2];

if($#ARGV < 2){ 
    print "There is not enough argument for BSSE calculation: input number of ghost atoms!\n";
    print "The execution is stopped.\n";
    exit(0);
}

$Bname="co2";

$fname=$fname_pre.".in";

open(IN,"<$fname");

$fA_a=$fname_pre."_a";
$fA_ab=$fname_pre."_ab";
$fB_b=$Bname."_b";
$fB_ab=$Bname."_ab";

# input file for parallel running
$inp1="c".$fA_a.".in";	
$inp2="b".$fA_ab.".in";	
$inp3="e".$fB_b.".in";	
$inp4="d".$fB_ab.".in";	
$inp5="a".$fname;	

system("mv $fname $inp5");

open(OUT0,">$inp1");
open(OUT1,">$inp2");
open(OUT2,">$inp3");
open(OUT3,">$inp4");

$i=0;
$tag=$i;
$flag="OFF";
$nskip=1;
@line_list=();

while($line=<IN>){
    @field=split(/\s+/,$line);
#    print "0:",$field[0],"1:",$field[1],"2:",$field[2],"\n";

    if($field[0] eq "\$molecule" and $flag eq "OFF"){
	$flag="ON";
	$tag=$i;
	$i++;
	print $line;
	print_file();
	next;
    }elsif($flag eq "ON"){
	if($i==$tag+1){
	    print $line;
	    print_file();
	}else{
	    $field[1]="@".$field[1];
	    $new_line="  ".$field[1]."\t".$field[2]."\t".$field[3]."\t".$field[4]."\n";
	    if($n<$Aatom){
		print OUT0 $line;
		print OUT1 $line;
		print OUT3 $new_line;
	    }elsif($n<$Aatom+$Batom){
		print OUT1 $new_line;
		print OUT2 $line;
		print OUT3 $line;
	    }else{
	    	print_file();    
	    }
	    print $line;
	    $n++;
	}
	$i++;
    }else {
    	print $line;
    	print_file();
    }
}

close(IN);
close(OUT0);	close(OUT1);	close(OUT2); 	close(OUT3);

sub print_file {
    print OUT0 $line;
    print OUT1 $line;
    print OUT2 $line;
    print OUT3 $line;
}
