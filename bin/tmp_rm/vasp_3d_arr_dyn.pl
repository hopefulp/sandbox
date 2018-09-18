#!/usr/bin/perl
## written by Joonho Park
## read POSCAR file 
## rearrangement of atoms for initial magnetic moment assign
## input: poscar_name atom_name a b c [direction priority]
#

use strict;
use warnings;
use Data::Dumper qw(Dumper);
use Data::Dump   qw(dump);
if($#ARGV < 3 or $ARGV[0] =~ /Usage/i){
    print "Error: input inputfile atom_name and lattice axis priority \n";
    print "Usage: $0 input_file Ni a b c \n";
    exit(1);
}

my $finp=shift @ARGV;
my $fout="CONTCAR";
print "INPUT OUTPUT fname: $finp & $fout\n";

my $atom_name=shift @ARGV;

open(IN,"<$finp");
open(OUT,">$fout");

##### entire range variable
my $d_limit=2.0;
#############################

my $i=0;
my $i_skip=0;	# skip atom before rearrangement
my $i_atom=0;	# count atom in rearrangement
my $n_skip=0;
my @atoms_arr=();
my @lattice_const=();
my $L_sort="N";
my ($j, $line, @line, $index, $i_species, @species_num, $natom_arrange );
my $dynamics="Y";
##### for output
my $debug="Y"; my $debug2="N";
my $output="Y";
if($output eq "Y"){ $debug="N";}

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
	if($i==2){ push (@lattice_const,$line[0]); print OUT $line,"\n"; print $line,"\n";next LINE;}
 	if($i==3){ push (@lattice_const,$line[1]); print OUT $line,"\n"; print $line,"\n";next LINE;}
	if($i==4){ push (@lattice_const,$line[2]); print OUT $line,"\n"; print $line,"\n";next LINE;}
	if($i==5){
	    print OUT $line."\n";
	    print  $line."\n";
	    $index=0;
	    ++$index until $line[$index] eq $atom_name or $index > $#line;
	    $i_species=$index;
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
	    $natom_arrange=$species_num[$i_species];
#	    print "@species_num"."\n";
	    next LINE;
	}
	if($line[0] !~ /^[\d\.]/){  
	    if($line[0] !~ /^sS/ and $dynamics eq "Y"){
		print OUT "Selective\n";
		print     "Selective\n";
	    }
	    print OUT $line."\n";
	    print  $line."\n";
	    next LINE; 
	}  # print $line[0]."\n";
	if($i_skip < $n_skip){ #cmp num_skip
	    if($dynamics eq "Y"){
		print OUT $line,"        T   T   T\n";
		print     $line,"        T   T   T\n";
	    }else{
	    	print OUT $line."\n";
	    	print $line."\n";
	    }
	    $i_skip++; next LINE;
	}
	if($i_atom<$natom_arrange){
	    push(@{$atoms_arr[$i_atom]}, @line);
	    $i_atom++;
	    next LINE;
	}
	#### get rearranged array ########
        &rearrange(\@ARGV,\@lattice_const,$natom_arrange,\@atoms_arr);
#	dump @atoms_arr;
	for($j=0;$j<$natom_arrange;$j++){
	    if($dynamics eq "Y"){
		print OUT "         ", join("   ",@{$atoms_arr[$j]}),"     T   T   T\n";
		print     "         ", join("   ",@{$atoms_arr[$j]}),"     T   T   T\n";
	    }else{
	    	print OUT "         ", join("   ",@{$atoms_arr[$j]}),"\n";
	    	print "         ",join("   ",@{$atoms_arr[$j]}),"\n";
	    }
	}
	$L_sort="Y";
	### add 1st line of new bloci
	if($dynamics eq "Y"){
            print OUT  $line,"    T   T   T\n";
            print      $line,"    T   T   T\n";
	}else{
            print OUT  $line;
            print      $line;
	}
    }else{
	if($dynamics eq "Y"){
	    chomp($line);
            print OUT $line,"     T   T   T\n";;
            print     $line,"     T   T   T\n";;
	}else{
            print OUT  $line;
            print      $line;
	}
    }
}    continue {
     $i++;
}
close(IN);
#print "skip $n_skip atoms and rearrange $natom_arrange\n";

#for($i=0;$i<$natom_arrange;$i++){
#    print "$atoms_arr[$i*3]  $atoms_arr[$i*3+1]  $atoms_arr[$i*3+2]\n";
#}
### Write with rearrangement
close(OUT);

##### SUBROUTINE ######
####  get rearranged atom list ####
sub rearrange { ### ref ref int ref
    my ($paxes, $plattice_cont, $natom, $coor_xyz)=@_;
    my @axes = @{$paxes};
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
    $i_axis[0]=get_i_axis($axes[0]);
    $i_axis[1]=get_i_axis($axes[1]); 
    $i_axis[2]=get_i_axis($axes[2]);

    ##### arrange following 1st priority in direction / make groups
    @group_up=get_sort($i_axis[0],$lattice_cons[$i_axis[0]],$natom, \@coord_xyz);
#    print Dumper \@coord_xyz;
    if($debug eq "Y") { print "1st group"; dump @group_up; }
    for($j=0,$if=0;$j<=$#group_up;$j++){
	if($j==0){ $in=$if;}
	else {	   $in=$if+1;}
	$if=$in+$group_up[$j];
#	print "divide groups from $in to $if \n";					#print -2
	@tmp_coord_arr=@coord_xyz[$in..$if];
#	dump @tmp_coord_arr;
	##### 2nd level sort in each group
	@group_mid=get_sort($i_axis[1],$lattice_cons[$i_axis[1]],$group_up[$j]+1,\@tmp_coord_arr);
#	dump @tmp_coord_arr;
	for($k=0,$ifc=0;$k<=$#group_mid;$k++){
            if($k==0){ $inc=$ifc;}
	    else {     $inc=$ifc+1;}
            $ifc=$inc+$group_mid[$k];
#            print "divide groups from $inc to $ifc \n";				# print -1
            @tmp_coord_arrc=@tmp_coord_arr[$inc..$ifc];
	    @group_low=get_sort($i_axis[2],$lattice_cons[$i_axis[2]],$group_mid[$k]+1,\@tmp_coord_arrc);
#	    dump @tmp_coord_arrc;							# print 0
	    splice(@tmp_coord_arr,$inc,$group_mid[$k]+1,@tmp_coord_arrc);
	}
	splice(@coord_xyz,$in,$group_up[$j]+1,@tmp_coord_arr);
    }
#    dump @coord_xyz;	
    @{$coor_xyz}=@coord_xyz;
    return 1;
}
##### devide into several groups
sub get_sort {
    my ($i_axis,$i_latt_const, $natom, $coord3)=@_;
    my (@sort_G2d, @sort_coord, @ngroup, @g_num);
    my ($j, $k, $inter, @ave, $L_insert, %sort, @tmp, @tmp2, @tmp_coord);
    my ($coord, $dist);
    @sort_G2d=(); @sort_coord=(); @ngroup=(); @g_num=();
    @ave=(); %sort=(); @tmp=(); @tmp2=(); @tmp_coord=();
#    print $coord3,"\n"; 
#    print Dumper ${$coord3}[0][0],"\n"; 
#    print "sort by $i_axis-axis with $natom atoms\n";			# print 4
    for my $i (0 .. $natom-1){
	$coord=${$coord3}[$i][$i_axis];
	#print $coord."\n";
	if($i==0){
	    push @{$sort_G2d[0]}, $i; 
	    $ave[0]=$coord;
	}else{
	    $L_insert="N";
	    #### loop for looking for the group
	    for($j=0;$j<=$#ave;$j++){
	        $dist=pbc(\$ave[$j], $coord, $i_latt_const);    #abs($ave[$j]-$coord);
	    	if($dist <= $d_limit){
		    push @{$sort_G2d[$j]}, $i;
		    $L_insert="Y";
		    last;
	    	}
	    }
	    if($L_insert eq "Y"){
		next;
	    }else{
		push @{$sort_G2d[$j]}, $i;
#	    	print "$i: sort_G2d [$j] [ $#{$sort_G2d[$j]} ] \n";
	 	$ave[$j]=$coord;
		next;
	    }
	}
    }
#    dump @sort_G2d;							# print 3
    #### making hash for group
    for($j=0;$j<=$#ave;$j++){
	$sort{$ave[$j]}= "@{$sort_G2d[$j]}";
    }
#    dump %sort;							# print 1
    #### sort group
    foreach my $x (sort {$a <=> $b} keys %sort){
	@tmp2=split(/\s+/,$sort{$x});
	push @g_num, $#tmp2;
	push @tmp, @tmp2;
    }
#    dump @tmp;								# print 2

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

