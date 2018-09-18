#!/usr/bin/perl
## written by Joonho Park
## part of "Vrun.sh"
## read POSCAR, metal  and modify
#
#

$poscar=$ARGV[0];
$metal=$ARGV[1];

open(IN,"< $poscar");
open(OUT,"> POSCAR");

LINE: while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($i==0){  print OUT "   $metal   $line[1]   $line[2]   $line[3]\n";
                print     "   $metal   $line[1]   $line[2]   $line[3]\n";
	next LINE;
    }
    if($i==5){
	if($line[0] !~ /[A-Z]/) {
            print "Error in $fin format: there are not atom symbols in $i+1 line\n";
            exit(0);
        }else{ 
	    print OUT "   $metal   $line[1]   $line[2]   $line[3]\n";
	    print     "   $metal   $line[1]   $line[2]   $line[3]\n";
 	}
	next LINE;
    }

    print OUT "$line\n";
    print     "$line\n";
} continue {
    $i++;
}

close (IN);
close (OUT);

