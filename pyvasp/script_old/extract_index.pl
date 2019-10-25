#!/usr/bin/perl
# written by Joonho Park
# read  POSCAR and extract atoms with index as xyz format
# $0 a.pos file.index

use lib '/qcfs/joonho/modules';
use Arith;
use Data::Dumper qw(Dumper);
use List::Util qw(sum);

if($#ARGV < 0){
    print "Usage:: $0 poscar file.index file.charge\n";
    print "        extract indices in mol format\n";
    print "        input POSCAR with Cartesian coordinates\n";
    exit(1);
}

$fposcar=$ARGV[0];
$f_index=$ARGV[1];
$f_chg=$ARGV[2];

$xyz_suf="xyz";
$mol_suf="mol";
$mol_ext_suf="ext_chg.mol";

@infile=split(/\./,$fposcar);
if($#infile == 1){
    $prefix=$infile[0];
}else{ print "Error(3):: there are many . in file name\n"; exit(3);}

# GET CO2 coordinates
open(IN1,"< $f_index");

### write distance
$outfile_d=$prefix.".dist";
open(OUTd,">$outfile_d");
@dist=();

### read INDEX file
### define number of output file from number of lines
$i=0;	# each indices in a line makes one out file
@katoms=();
@n_katoms=();
@extract_katom=(M, O, C, O);
@lattice_a=();
@lattice_b=();
@lattice_c=();
while($index=<IN1>){
    chomp($index);
    @index=split(/\s+/,$index);
    if($index[0] eq ""){shift(@index);}	

    $outfile=$prefix."-$i".".$mol_suf";
    $outfile2=$prefix."-$i".$mol_ext_suf;
    open(OUT,">$outfile");
    open(OUT2,">$outfile2");

    ### initial calculation
    if($i==0){
	$n_extract_atom=$#index+1;
    }

     #print "$n_extract_atom\n $outfile\n";
     #print OUT "$n_extract_atom\n $outfile\n";
     print "\n\$molecule\n\t0\t1\n";
     print OUT "\n\$molecule\n\t0\t1\n";
     print OUT2 "\n\$molecule\n\t0\t1\n";
   
    ### get charge for each index line
    open(IN3,"< $f_chg");
    $i_chg=0;
    while($line=<IN3>){
	if($i_chg<$i) { next; }
	elsif($i_chg == $i){
            chomp($line);
            @line=split(/\s+/,$line);
            if($line[0] eq ""){shift(@line);}
	    $chg=$line[0];
	    last;
	}else{ print "Error{5):: Algorithm error in getting charge\n"; exit(5);}
    } continue {
    	$i_chg++;
    }
    close(IN3);

    ### read POSCAR through 1 time for each index line
    $j=0;
    $tag="NO";
    $natom=0;
    @coord=();

    open(IN2,"< $fposcar");
    while($line=<IN2>){     # in the order of M O O C
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
	
        #print join(" ",@n_katoms),"\n";
 	$a=$j-8;	# a is atom index
	# if the atom index in poscar is the same as extract index
	for($k=0;$k<$n_extract_atom;$k++){
	    if($a==$index[$k]){
		@{$coord[$natom]}=($line[0], $line[1], $line[2]);
		#print "here: ",join(" ",@{$coord[$natom]}),"\n";
		$natom++;
		last;
	    }
	}
	#print "$natom ? $n_extract_atom\n";
	if($natom==$n_extract_atom){
	    #print Dumper \@coord;
	    if($index[1] < $index[3]){  # IN the order of M O1 C O2
		@tmp=@{$coord[2]};	# if		  M O1 O2 C
		@{$coord[2]}=@{$coord[3]};
		@{$coord[3]}=@tmp;
	    }else{			# if		  M O2 O1 C
		@tmp=@{$coord[2]};
		@{$coord[2]}=@{$coord[3]};
		@{$coord[3]}=@{$coord[1]};
		@{$coord[1]}=@tmp;
	    }
	    #print Dumper \@coord;
	    #$k=0;
	    #print "$katoms[$k] ",join("  ",@{$coord[$k]}),"\n";
	    for($k=1;$k<$n_extract_atom;$k++){
		#here distance test
		#dist_lattice_z: move atom to near location
		@dist_check=&Arith::dist_lattice_z(\@lattice_a,\@lattice_b,\@lattice_c,$coord[$k-1],$coord[$k]);
		#print "tag: $dist_check[0] $dist_check[1]\n";
		print "$extract_katom[$k] ",join("  ",@{$coord[$k]}),"\n";
		print OUT "$extract_katom[$k] ",join("  ",@{$coord[$k]}),"\n";
		print OUT2 "$extract_katom[$k] ",join("  ",@{$coord[$k]}),"\n";
		if($k==1){
		    push @dist, $dist_check[1];
		    print OUTd "$dist_check[1]\n";
		}
	    }

	    print "\$end\n";
	    print OUT "\$end\n\n";
	    print OUT2 "\$end\n\n";

	    ##### distanct test for periodic cell

	    last;
	}
	next;
    } continue {
    	$j++;
    }

    print OUT2 "\n\$external_charges\n",join("  ",@{$coord[0]}), "   $chg\n\$end\n\n";

    close(IN2); 
    close(OUT);
    close(OUT2);

} continue {
    $i++;
}

print "dist mean= ", mean(@dist),"\n";
print OUTd "dist mean= ", mean(@dist),"\n";

close(OUTd);
close(IN1); 

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

