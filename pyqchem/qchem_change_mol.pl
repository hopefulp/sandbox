#!/usr/bin/perl
# written by Joonho Park

$fname=$ARGV[0];
$jobname=$ARGV[1];
if($#ARGV < 1){
    print "less argument: can't change filename\n";
    exit(0);
}

@fname=split(/\./,$fname);
$newfile=$fname[0]."_$jobname.inp";
#$newfile=$fname[0]."_$jobname.$fname[1]";

open(IN,"<$fname");
open(OUT,">$newfile");
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
system("qchem $inp1 $out1 &");
system("qchem $inp2 $out2 &");
system("qchem $inp3 $out3 &");
system("qchem $inp4 $out4 &");
system("qchem $inp5 $out5 &");

sub print_file {
    print OUT0 $line;
    print OUT1 $line;
    print OUT2 $line;
    print OUT3 $line;
}
