#!/usr/bin/perl
# written by Joonho Park

if($#ARGV < 2){ 
    print "There is not enough argument for BSSE calculation: input number of ghost atoms!\n";
    print "The execution is stopped.\n";
    exit(0);
}

$fname_pre=shift(@ARGV);
$Aatom=shift(@ARGV);
$Batom=shift(@ARGV);

$Bname="co2";

$fname=$fname_pre.".in";

open(IN,"<$fname");

$fA_a=$fname_pre."_a";
$fA_ab=$fname_pre."_ab";
$fB_b=$Bname."_b";
$fB_ab=$Bname."_ab";

$inp1="c".$fA_a.".in";	$out1="c".$fA_a.".out";
$inp2="b".$fA_ab.".in";	$out2="b".$fA_ab.".out";
$inp3="e".$fB_b.".in";	$out3="e".$fB_b.".out";
$inp4="d".$fB_ab.".in";	$out4="d".$fB_ab.".out";
$inp5="a".$fname;	$out5="a".$fname_pre.".out";

open(OUT0,">$inp1");
open(OUT1,">$inp2");
open(OUT2,">$inp3");
open(OUT3,">$inp4");
open(OUT4,">$inp5");

$i=0;
$tag=$i;
$flag="OFF";
$nskip=1;
@line_list=();

$new_ch=0;
$new_mul=1;

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
	    print OUT0 $line;
	    print OUT1 $line;
	    print OUT2 "\t$new_ch\t$new_mul\n";
	    print OUT3 "\t$new_ch\t$new_mul\n";
	    #print_file();
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
system("mv $fname $inp5");

# when sp of complex has been done
system("rm $inp5");

# for serial job
#system("qchem $inp1 $out1 &");
#system("qchem $inp2 $out2 &");
#system("qchem $inp3 $out3 &");
#system("qchem $inp4 $out4 &");
#system("qchem $inp5 $out5 &");

sub print_file {
    print OUT0 $line;
    print OUT1 $line;
    print OUT2 $line;
    print OUT3 $line;
}
