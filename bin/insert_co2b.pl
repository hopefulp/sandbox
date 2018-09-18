#!/usr/bin/perl
# written by Joonho Park
# read  bare mof POSCAR	 and add CO2 extracted from mof-co2 POSCAR and output file


if($#ARGV < 1){
    print "Error: input TWO POSCAR's ";
    exit(0);
}

$fmof=$ARGV[0];
$fco2=$ARGV[1];

if($#ARGV >= 2){
    $newfile=$ARGV[2];
}else{
    $newfile=$fmof.".new";
}

print "New file is $newfile\n";

# GET CO2 coordinates
open(INc,"< $fco2");

$line_atom_number="    6    30    30    6\n";
$nOplus=12;
$nCplus=6;
$nOstart=24;	# 6(Zn) + 18(O)
$nCstart=60;	# 6(Zn) + 30(O) + 24(C)
$tag="OFF";
$natom_total=72;
$selective_dynamics="NO";
$i=0;
$natom=0;
### read three atoms's coordinates
while($line=<INc>){
#    if($i<5){	$i++;	next; }
    @field=split(/\s+/,$line);
    if($field[0] eq ""){ shift(@field);}
#    print @field[0]."\n";
    if($tag eq "OFF" and $field[0] ne "Direct"){
	$i++;    	next;
    }elsif($tag eq "OFF" and $field[0] eq "Direct"){
	$tag="ON";	$i++;	next;
    }
#    print "$i= ", $tag,$natom,"\n";
    if($natom < $nOstart){ $natom++; next;}
    elsif($natom < $nOstart+$nOplus){ 
	$coord[$natom][0]=$field[0];
	$coord[$natom][1]=$field[1]; 
	$coord[$natom][2]=$field[2];
#	print $coord[$natom][0];
	$natom++;	next;
    }
    
    if($natom < $nCstart){$natom++; next;}
    elsif($natom < $nCstart+$nCplus){
 	$coord[$natom][0]=$field[0];
        $coord[$natom][1]=$field[1];
        $coord[$natom][2]=$field[2];
#	print $coord[$natom][0];
       $natom++;       next;
    }
    
}
close(INc); 
#print 
#if($natom != $nOplus+$nCplus){
#    print "not match the number of CO2\n";
#    exit(0);
#}

$i=0;
$natom=0;
$tag="OFF";
open(INm,"< $fmof");
open(OUT,"> $newfile");
while($line=<INm>){
#    if($i<5){  $i++;   next; }
    @field=split(/\s+/,$line);
    if($field[0] eq ""){ shift(@field);}
    print $field[0],"\n";
    if($tag eq "OFF"){
        if    ($field[0] eq "Direct"){
            print OUT $line;
	    $tag="ON";      $i++; 
	}elsif($field[0] eq "Selective"){
            print OUT $line;
	    $selective_dynamics="YES"; $i++;
	}elsif($field[0] == 6){
	    print OUT $line_atom_number;
	    $i++; 
        }else{  print OUT $line;}
	next;
    }

    if($natom < $nOstart){  print OUT $line; $natom++; next;}
    elsif($natom < $nOstart+$nOplus){
	for($j=0;$j<$nOplus;$j++){
            print OUT "  ",$coord[$natom][0];
            print OUT "  ",$coord[$natom][1];
            print OUT "  ",$coord[$natom][2];
	    if($selective_dynamics eq "YES"){ 	print OUT "   T   T   T\n";}
	    else{ 					print OUT "\n";}
	    print "Write new atom: ",$coord[$natom][0]."\n";
	    $natom++;
	} 
    	print $natom,"check\n";
	print OUT $line; $natom++; next; # write the present atom after adding new atoms
    }
    if($natom < $nCstart){  print OUT $line; $natom++; next;}
    elsif($natom < $nCstart+$nCplus){
        for($j=0;$j<$nCplus;$j++){
            print OUT "  ",$coord[$natom][0];
            print OUT "  ",$coord[$natom][1];
            print OUT "  ",$coord[$natom][2];
            if($selective_dynamics eq "YES"){   print OUT "   T   T   T\n";}
            else{                                       print OUT "\n";}
	    print "Write new atom: ",$coord[$natom][0]."\n";
	    $natom++;
        }
    	print $natom,"check\n";
	print OUT $line; $natom++; next; # write the present atom after adding new atoms
     }

    if($natom < $natom_total) {print OUT $line; $natom++;}
}
close(INm);
close(OUT); 
