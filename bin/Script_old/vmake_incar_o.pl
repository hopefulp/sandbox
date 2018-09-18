#!/usr/bin/perl
## written by Joonho Park
## Usage:: vmake_incar incar.FM
#
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
open(OUT,"> INCAR");

LINE: while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($line[0] eq "MAGMOM"){
	printf OUT ("MAGMOM = %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f 999*0\n", $magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal}) ;
	printf     ("MAGMOM = %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f 999*0\n", $magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal},-$magmom{$metal}) ;
    next LINE;
    }

    print OUT "$line\n";
    print     "$line\n";
}   continue {
    $i++;
}

close(IN);
close(OUT);


