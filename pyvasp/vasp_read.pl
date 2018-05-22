#!/usr/bin/perl
## written by Joonho Park
## read POSCAR file 
## rearrangement of atoms for initial magnetic moment assign
## input: atom_name a b c [direction priority]
#

if($#ARGV < 3 or $ARGV[0] =~ /Usage/i){
    print "Error: input atom_name and lattice axis priority \n";
    print "Usage: $0 Ni a b c \n";
    exit(1);
}

$finp="POSCAR";
$fout="CONTCAR";
print "INPUT OUTPUT fname: $finp & $fout\n";

$atom_name=shift @ARGV;

open(IN,"<$finp");
open(OUT,">$fout");


$i=0;
$i_skip=0;	# skip atom before rearrangement
$i_atom=0;	# count atom in rearrangement
$n_skip=0;
@atoms_arr=();
LINE: while($line=<IN>){
	next LINE if($i<5);
	
	chomp($line);
	@line=split(/\s+/,$line);
	if($line[0] eq ""){ shift(@line);}
	if($i==5){
	    $index=0;
	    ++$index until $line[$index] eq $atom_name or $index > $#line;
#	    $index = grep { $line[$_] eq $atom_name } 0..$#line; # not working
	    $i_species=$index;
	    print "$line\n";
	    print "Searched Atom $atom_name index: $index\n";
	    next LINE;
	}
	if($i==6){
	    @species_num=@line;
	    for($j=0;$j<$i_species;$j++){
		$n_skip+=$species_num[$j];
	    }
	    $natom_arrange=$species_num[$i_species];
	    print "@species_num"."\n";
	    next LINE;
	}
	if($line[0] !~ /^[\d\.]/){  next LINE; }  # print $line[0]."\n";
	if($i_skip < $n_skip){ #cmp num_skip
	    $i_skip++; next LINE;
	}
	if($i_atom<$natom_arrange){
	    push(@atoms_arr, @line);
	    $i_atom++;
	    next LINE;
	}
}    continue {
     $i++;
}
close(IN);
#print "skip $n_skip atoms and rearrange $natom_arrange\n";

for($i=0;$i<$natom_arrange;$i++){
    print "$atoms_arr[$i*3]  $atoms_arr[$i*3+1]  $atoms_arr[$i*3+2]\n";
}





close(OUT);
