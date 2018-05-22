#!/usr/bin/perl
# written by Joonho Park
# read POSCAR and extract atoms as index using anchor
# $0 a.pos file.index

use lib '/qcfs/joonho/modules';
use Arith;
use Data::Dumper qw(Dumper);
use List::Util qw(sum);

if($#ARGV < 0){
    print "Usage:: $0 poscar \" indices \" \n";
    print "        extract indices in the order of connection\n";
    print "        input POSCAR with Cartesian coordinates\n";
    exit(1);
}

$fposcar=shift(@ARGV);
$anchor=shift(@ARGV);
$atom_string=shift(@ARGV);

##### read poscar
@infile=split(/\./,$fposcar);
if($#infile == 1){
    $prefix=$infile[0];
}else{ print "Error(3):: there are many . in file name\n"; exit(3);}

#### read pivot atoms
@anchor=split(/\s+/,$anchor);
print join(" ",@anchor),"\n";

#### extract atoms
@atom_string=split(/\s+/,$atom_string);
print join(" ",@atom_string),"\n";

### write distance
$outfile_d=$prefix.".dist";
$outfile=$prefix.".index_dist";
open(OUT,">$outfile");
@dist=();

@katoms=();
@n_katoms=();
@lattice_a=();
@lattice_b=();
@lattice_c=();
@atoms=();

$low_limit1=2.2;
$high_limit1=2.9;
$low_limit2=1.1;
$high_limit2=1.6;

open(IN,"< $fposcar");
### read POSCAR and save atoms coordinates
$j=0;
while($line=<IN>){
	if($j<2){ next; }
	chomp($line);
    	@line=split(/\s+/,$line);
    	if($line[0] eq ""){shift(@line);}
	
	if($j==2) {@lattice_a=@line; next;}
	if($j==3) {@lattice_b=@line; next;}
	if($j==4) {@lattice_c=@line; next;}
	if($j==5) {@katoms=@line;   next;}
	if($j==6) {@n_katoms=@line; next;}
	if($j==7){ 
    	    if($line[0] =~ /^[dD]/) {print "Error(2)::Input file is supposed to have cart coordinates\n"; exit;
	    }else{ next;}
	}
 	$a=$j-8;	# a is atom index
	# save atoms in @atoms
	#push @atoms, @line;	# as 1d array
	for($i=0;$i<3;$i++){
	    $atoms[$a][$i]=$line[$i];
	}
}continue {
    	$j++;
}
close(IN);

$tot_atom=sum(@n_katoms);
#print "$tot_atom\n";

if($tot_atom != $a+1){
    print "Read poscar error\n";
    exit(5);
}
#print Dumper \@atoms;

@pivot=();
@conn=();
#### find connected atoms in @atoms
#### loop for anchor atoms in command line argument
for($i=0;$i<=$#anchor;$i++){
    # save pivot atom_kind position in poscar
    @pivot=@{$atoms[$anchor[$i]]};
    @index=();
    print Dumper \@pivot;
    #### loop for atom_string in command line argument
    for($j=1;$j<=$#atom_string;$j++){
	if($j==1){$low_limit=$low_limit1; $high_limit=$high_limit1;}
	else	{ $low_limit=$low_limit2; $high_limit=$high_limit2;}
	$conn_atom=$atom_string[$j];
	$k_find_tag="NO";
	### find the connected atom species
	for($k=0;$k<=$#katoms;$k++){
	    if($conn_atom eq $katoms[$k]){
		$i_kind=$k;
		$k_find_tag="YES";
		last;
	    }
	}
	print "species = $i_kind\n";
	if($k_find_tag eq "NO"){ print "Error in finding atom species\n"; exit(11);}
	### count atoms pre-resided
	$pre_atoms=0;
	for($k=0;$k<$i_kind;$k++){
	    $pre_atoms+=$n_katoms[$k];
	}
	print " number of pre-atoms: $pre_atoms\n";
 	### find the connected atoms
 	LOOP2: for($k=0;$k<$n_katoms[$i_kind];$k++){
	    $n=$pre_atoms+$k;
	    for($l=0;$l<=$#index;$l++){
		if($index[$l] == $n){ next LOOP2; } 
	    }
	    #if($j==2) {print "$n atom coddord: @{$atoms[$n]}\n";}
	    @dist_check=&Arith::dist_lattice(\@lattice_a,\@lattice_b,\@lattice_c,\@pivot,$atoms[$n]);
	    #if($j==2) {print "$n: $dist_check[0] $dist_check[1]\n";}
	    if($low_limit<$dist_check[1] and $dist_check[1]<$high_limit){
	  	$index[$j-1]=$n;
		@pivot=@{$atoms[$n]};
		print "index :$index[$j-1]  $dist_check[1]\n";
		if($j==1) { push @dist, $dist_check[1];}
		last LOOP2;
	    }
	}
    }
    print OUT "$anchor[$i] ",join(" ",@index), "  $dist[$i] \n";
}


print "dist mean= ", mean(@dist),"\n";
print OUT "dist mean= ", mean(@dist),"\n";

close(OUT);

###################### Subroutine
sub mean { return @_ ? sum(@_)/@_ : 0 }

sub get_atom {
    my ($n, $k_atoms, $n_k_atoms)=@_;
    my ($i, $nsum_atom, @atom_species, @n_atom_species);
    @atom_species=@{$k_atoms};
    @n_atom_species=@{$n_k_atoms};

    #confirm referece and pointer
    #print join(" ",@atom_species),"\n";
  
    $nsum_atom=0; 
    for($i=0;$i<$#atom_species;$i++){
	$nsum_atom+=$n_atom_species[$i];
	if($n <= $nsum_atom){ return $atom_species[$i]; }
	else {next;}
    }
    print "Error(10):: algorithm error in sub atom\n"; exit(10);
}

