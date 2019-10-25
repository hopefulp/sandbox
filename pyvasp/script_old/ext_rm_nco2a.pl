#!/usr/bin/perl
# written by Joonho Park
# read  mof-co2 POSCAR 
# with direct coordinate ./a.pl poscar D
# with cartesian coordinate and extract co2 in .6co2.pos and .6co2.xyz; ./a.pl poscar (C)

if($#ARGV < 0){
    print "Error: it needs input file ";
    exit(0);
}
$fin=$ARGV[0];

@rm_atom=(25,26,27,28,29,30,33,34,35,36,61,62,63,65,66);
$nmol=5;
$new_suffix="ads1co2";
$fout=$fin.".$new_suffix";

# GET New MOF with 1CO2 coordinates
open(IN,"< $fin");
open(OUT,"> $fout");

@extr_atom=qw( O C );
@natom_ea=qw( 2 1 );

$test="OFF";	# atom copy test
$tag="OFF"; 	# atom copy tag
$ikind=0;

@out_line=();
$ntatoms=0;
$i=0;
$iatom=0;
### read three atoms's coordinates
LINE: while($line=<IN>){
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}

    if($i==0){  print OUT join("  ",@extr_atom), "\n";
		print     join("  ",@extr_atom), "\n";
		push(@out_line,join("  ",@extr_atom)."\n"); next LINE;}
    if($i==1){  print OUT $line."\n"; 
	       	print 	  $line."\n";
		push(@out_line,$line."\n"); next LINE;}
    if($i<5){   print OUT $line."\n"; 
		print 	  $line."\n"; 
		push(@out_line,$line."\n"); next LINE;}
    if($i==5){
        if($line[0] !~ /[A-Z]/) {
            print "Error in $fin format: there are not atom symbols in $i+1 line\n";
            exit(0);
        }else{ for($j=0;$j<=$#line;$j++){push(@atoms,$line[$j]);} }
	print OUT join("  ",@atoms)."\n";
	print     join("  ",@atoms)."\n";
        next LINE;
    }
    if($i==6){
        for($j=0;$j<=$#line;$j++){
            $ntatoms+=$line[$j];
	    $iextr=-1;
	    for($k=0;$k<=$#extr_atom;$k++){
		if($atoms[$j] eq $extr_atom[$k]){
		    $iextr=$k;
		    last;
	     	}
  	    }
	    if($iextr != -1){
		$n_extr_atom=$natom_ea[$iextr]*$nmol;
	    }else{ $n_extr_atom=0;}
            $natom_kind[$j]=$line[$j]-$n_extr_atom;	# number of atoms of each kind
        }
	print OUT join("  ",@natom_kind),"\n";
	print     join("  ",@natom_kind),"\n";
	push(@out_line,join("  ",@natom_extracted)."\n");
        next LINE;
    }
    next LINE if($line[0] =~ /^[sS]/);	# for selected dynamics
    if($line[0] =~ /^[dD]/){
   	print OUT $line."\n"; 
   	print     $line."\n"; 
	next LINE;
    }

    # AFTER knowing 
    # WRITE or SKIP test 
    $check_tag="N";
    for($j=0;$j<=$#rm_atom;$j++){
	if($iatom+1 == $rm_atom[$j]){
	    $check_tag="Y";
	    last;
	}
    }
    if($check_tag eq "Y"){
	$iatom++;
	next;
    }else{
	print OUT $line."\n";   
	print     $line."\n";   
	$iatom++;
    }

} continue {
    $i++;
}


close(IN); close(OUT);
