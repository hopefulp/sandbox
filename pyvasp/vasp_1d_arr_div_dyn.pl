#!/usr/bin/perl
## written by Joonho Park
## read POSCAR file 
## rearrangement of atoms for initial magnetic moment assign
## input: poscar_name atom_name a averages[groups]
#

use strict;
#use warnings;
use Data::Dumper qw(Dumper);
use Data::Dump   qw(dump);
if($#ARGV < 3 or $ARGV[0] =~ /Usage/i){
    print "Error: input inputfile atom_name and lattice axis priority \n";
    print "Usage: $0 inputfile dynamic_option Ni a[bc] ave1 ave2 ave3 ... \n";
    exit(1);
}

my $finp=shift @ARGV;
my $fout=$finp.".1";
print "INPUT OUTPUT fname: $finp & $fout\n";

my $dynamics=shift @ARGV;
my $atom_name=shift @ARGV;
my $axis=shift @ARGV;
open(IN,"<$finp");
open(OUT,">$fout");

die "dynamics option should be \"N\" or \"Y\"" if($dynamics ne "Y" and $dynamics ne "N");

##### entire range variable
my @d_sect=@ARGV;
#############################

my $i=0;
my $i_skip=0;	# skip atom before rearrangement
my $i_atom=0;	# count atom in rearrangement
my $n_skip=0;   # number of atom species to be skipped
my @atoms_2b_arr=();
my @lattice_const=();
my $L_sort="N";
my ($j, $line, @line, $index, $i_species, @species_num, $natom_2b_arr, $pos_atom, $last_atom );
$last_atom="N";

##### for output
my $output="Y";
my $debug="N"; my $debug1="N"; my $debug_2dm="N"; my $debug3="N";
if($output eq "N"){
 $debug="Y";  $debug_2dm="Y";
 $debug1="Y"; $debug3="N";
}

LINE: while($line=<IN>){
    if($L_sort eq "N"){
	if($i<2){
	    print OUT $line;
	    print  $line;
	    next LINE;
	}
	chomp($line);
	@line=split(/\s+/,$line);
	if($line[0] eq ""){ shift(@line);}
	#print "* $line[0] \t  **$line[1] \t  ***$line[2] \t  ****$line[3]\n";


	if($i==2){ push (@lattice_const,$line[0]); print OUT $line,"\n"; print $line,"\n";next LINE;}
 	if($i==3){ push (@lattice_const,$line[1]); print OUT $line,"\n"; print $line,"\n";next LINE;}
	if($i==4){ push (@lattice_const,$line[2]); print OUT $line,"\n"; print $line,"\n";next LINE;}
	if($i==5){
	    print OUT $line."\n";
	    print  $line."\n";
	    $index=0;
	    ++$index until $line[$index] eq $atom_name or $index > $#line;
	    $i_species=$index;
	    $pos_atom=$line[$i_species];
#	    print "$line\n";
#	    print "Searched Atom $atom_name index: $index\n";
	    next LINE;
	}
	if($i==6){
	    print OUT $line."\n";
	    print  $line."\n";
	    @species_num=@line;
	    for($j=0;$j<$i_species;$j++){
		$n_skip+=$species_num[$j];
	    }
	    $natom_2b_arr=$species_num[$i_species];
	    if($debug eq "Y") {print "POSCAR in: $pos_atom $species_num[$i_species] rearrange\n";}
#	    print "@species_num"."\n";
	    next LINE;
	}
	if($line[0] !~ /^[\d\.-]/){  
	    if($line[0] !~ /^sS/ and $dynamics eq "Y"){
		print OUT "Selective\n";
		print     "Selective\n";
	    }
	    print OUT $line."\n";
	    print  $line."\n";
	    next LINE; 
	}   
	if($debug1 eq "Y"){ print $line[0] ,"\n";}
	if($i_skip < $n_skip){ #cmp num_skip
	    if($dynamics eq "Y"){
		print OUT $line,"        T T T\n";
		if($output eq "Y") {print     $line,"        T T T $i_skip\n";}
	    }else{
	    	print OUT $line."\n";
	    	if($output eq "Y") {print $line." $i_skip \n";}
	    }
	    $i_skip++; 
	    if($i_skip==$n_skip and $debug eq "Y"){ print "natom before arr: $i_skip, 2b_arr: $natom_2b_arr\n";}
	    next LINE;
	}
	if($i_atom<$natom_2b_arr){
	    push(@{$atoms_2b_arr[$i_atom]}, @line);
	    $i_atom++;
	    if($i_atom != $natom_2b_arr) {next LINE;}
	}
	#### get rearranged array ########
	if($debug1 eq "Y") {print "start rearrange () \n"; }
        &rearrange(\$axis,\@lattice_const,$natom_2b_arr,\@atoms_2b_arr);
#	dump @atoms_2b_arr;
	for($j=0;$j<$natom_2b_arr;$j++){
	    if($dynamics eq "Y"){
		print OUT "         ", join("   ",@{$atoms_2b_arr[$j]}),"     T T T\n";
		if($output eq "Y") {print     "         ", join("   ",@{$atoms_2b_arr[$j]}),"     T T T $pos_atom $j\n";}
	    }else{
	    	print OUT "         ", join("   ",@{$atoms_2b_arr[$j]}),"\n";
	    	if($output eq "Y") {print "         ",join("   ",@{$atoms_2b_arr[$j]}),"  $pos_atom \n";}
	    }
	}
	$L_sort="Y";
    }else{
	if($dynamics eq "Y"){
	    chomp($line);
            print OUT $line,"     T T T\n";;
            print     $line,"     T T T\n";;
	}else{
            print OUT  $line;
            print      $line;
	}
    }
}    continue {
     $i++;
}
print "number of line of inputfile: $i\n";

close(IN);
#print "skip $n_skip atoms and rearrange $natom_2b_arr\n";

#for($i=0;$i<$natom_2b_arr;$i++){
#    print "$atoms_2b_arr[$i*3]  $atoms_2b_arr[$i*3+1]  $atoms_2b_arr[$i*3+2]\n";
#}
### Write with rearrangement
close(OUT);

##### SUBROUTINE ##########################
####  rearrange the selected atom list ####
sub rearrange { ### ref ref int ref
    my ($paxes, $plattice_cont, $natom, $coor_xyz)=@_;
    my $paxis = ${$paxes};
    my @lattice_cons = @{$plattice_cont};
    my @coord_xyz = @{$coor_xyz};
    # for index sort
    my (    @group_up,  @i_sort, @i_axis );
    my ($j, @group_mid, @tmp_coord_arr,  $in,  $if);
    my ($k, @group_low, @tmp_coord_arrc, $inc, $ifc);
    @group_up=();  @i_sort=(); @i_axis=();
    @group_mid=(); @tmp_coord_arr=();
    @group_low=(); @tmp_coord_arrc=();
    ### index 1, 2, 3 has priority in itself
    ##### obtain axis priority in number
#    $i_axis[0]=&get_i_axis($axes[0]);
#    $i_axis[1]=&get_i_axis($axes[1]); 
#    $i_axis[2]=&get_i_axis($axes[2]);
    $i_axis[0]=&get_i_axis($paxis);

    ##### Arrange atoms following 1st priority in direction / make groups
#    dump @coord_xyz; exit;
    @group_up=&get_sort($i_axis[0],$lattice_cons[$i_axis[0]],$natom, \@coord_xyz);
#    print Dumper \@coord_xyz;
    if($debug eq "Y") { print "1st group: ", join("  ", @group_up),"\n"; }
    @{$coor_xyz}=@coord_xyz;
    return 1;
}
##### devide into several groups and return the rearranged list
sub get_sort {
    my ($i_axis,$i_latt_const, $natom, $coord3)=@_;
    my (@sort_G2d, @sort_coord, @ngroup, @g_num);
    my ($j, $k, $inter, $L_insert, %sort, @tmp, @tmp2, @tmp_coord);
    my ($coord, $dist);
    # Group 1-dimension 
    @sort_G2d=(); @sort_coord=(); @ngroup=(); @g_num=();
    %sort=(); @tmp=(); @tmp2=(); @tmp_coord=();
#    print $coord3,"\n"; 
#    print Dumper ${$coord3}[0][0],"\n"; 
#    print "sort by $i_axis-axis with $natom atoms\n";			# print 4
    if($debug1 eq "Y") {print "start get_sort() \n";}
    for my $i (0 .. $natom-1){
	$coord=${$coord3}[$i][$i_axis];
	if($debug1 eq "Y"){print "$coord  $d_sect[0]\n";}
	$L_insert="N";
	for($j=0;$j<=$#d_sect;$j++){
	    if($coord <= $d_sect[$j]){
		push @{$sort_G2d[$j]}, $i;
		$L_insert="Y";
		last;
	    }	
	    if($L_insert eq "N" and $j==$#d_sect){
	    	push @{$sort_G2d[$j+1]}, $i;
	    }
	}
    }
    if($debug eq "Y"){
	print "in &get_sort: sort_G2d \n";
    	dump @sort_G2d;
    }
    # transform 2d-array of @sort_G2d into 1d of @tmp
    for($j=0;$j<=$#d_sect+1;$j++){
	@tmp2=@{$sort_G2d[$j]};
	push @g_num, $#tmp2+1;
	push @tmp, @tmp2;
    }
    if($debug eq "Y"){ dump @tmp;}
    # change atom order following @tmp
    for($j=0;$j<=$#tmp;$j++){
        push @tmp_coord, @{$coord3}[$tmp[$j]];
    }
    @{$coord3}=@tmp_coord;
#    print Dumper @{$coord3};
    return @g_num;
}

sub pbc {
    my ($pivot, $coord, $lattice_const)=@_; # pivot(ave) is reference to be changed
    my ($dist, $dist2);

    $dist=abs($$pivot-$coord);
    $dist2=abs($lattice_const-$dist);
    if($dist2 < 2.0){ 
	$dist=$dist2;
	#### to move edge atoms to first group, $$pivot should be around 0
	if($coord < $$pivot) { $$pivot=$coord; }
#	print "dist2= $dist2: $lattice_const - $dist \n"; 
    }
    return $dist;
}

sub get_i_axis {
    my ($axis)=@_;

    if    ($axis eq "a"){ return 0;}
    elsif ($axis eq "b"){ return 1;}
    elsif ($axis eq "c"){ return 2;}
    else  {print "Error in axis priority\n";}
    exit(2);
}

