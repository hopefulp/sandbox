#!/usr/bin/perl
## written by Joonho Park
## Usage:: vmake_incar incar.FM
#	   this change magnetic moment value with the same sign
#

$incar=$ARGV[0];
$metal=$ARGV[1];

%magmom=(
'Sc'=>2 ,
'Ti'=>3 ,
'V'=>4.5 ,
'Cr'=>6 ,
'Mn'=>7.5 ,
'Fe'=>6 ,
'Co'=>4.5 ,
'Ni'=>3 ,
'Cu'=>2
);

open(IN,"< $incar");
open(OUT,"> INCAR.mag");

LINE: while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($line[0] eq "MAGMOM"){
	shift(@line);
	if($line[0] = "="){
	    shift(@line);
	}
	print OUT "MAGMOM = ";
	print     "MAGMOM = ";
	for($m=0;$m<$#line;$m++){
	    if($line[$m] lt 0){
		$magmom=-$magmom{$metal};
	    }else{
		$magmom=$magmom{$metal};
	    }
	    print OUT "$magmom  ";
	    print     "$magmom  ";
	}
	print OUT "$line[$#line]\n";
	print     "$line[$#line]\n";
    next LINE;
    }

    print OUT "$line\n";
    print     "$line\n";
}   continue {
    $i++;
}

close(IN);
close(OUT);


