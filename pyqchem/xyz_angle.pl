#!/usr/bin/perl
use Math::Trig;

$fin=$ARGV[0];

@CO2_O= qw( 25 26 27 28 29 30 31 32 33 34 35 36);
@CO2_C= qw( 61 62 63 64 65 66);

open(IN,"<$fin");

@O_coord=();
@C_coord=();
$i=0;
$index_atom=0;
while($line=<IN>){
    next if($i<2);

    chomp($line);
    @line=split(/\s+/, $line);
    if($line[0] eq ""){ shift(@line);}
    if($i==0 and $line[0] =~ /\D/){
	print "Error: In xyz format, first line is expected to be a number!\n";
	next;
    }
    
    $index_atom++;
#    print $index_atom,"\n"; 
    if(25 <= $index_atom and $index_atom <= 36){
	for($j=0;$j<3;$j++){
	    push(@O_coord,$line[$j+1]);
	}
	next;
    }elsif(61 <= $index_atom and $index_atom <= 66){
	for($j=0;$j<3;$j++){
            push(@C_coord,$line[$j+1]);
        }
        next;
    }
} continue {
    $i++; 
}

close(IN);

$sum_angle=0;
for($i=0,$j=0;$i<=$#O_coord;$i+=6,$j+=3){
    $dx1=$O_coord[$i+0]-$C_coord[$j+0];
    $dy1=$O_coord[$i+1]-$C_coord[$j+1];
    $dz1=$O_coord[$i+2]-$C_coord[$j+2];
    $dx2=$O_coord[$i+3]-$C_coord[$j+0];
    $dy2=$O_coord[$i+4]-$C_coord[$j+1];
    $dz2=$O_coord[$i+5]-$C_coord[$j+2];
    
    $dist1=sqrt($dx1*$dx1+$dy1*$dy1+$dz1*$dz1);
    $dist2=sqrt($dx2*$dx2+$dy2*$dy2+$dz2*$dz2);
    $inner_pro=$dx1*$dx2+$dy1*$dy2+$dz1*$dz2;

    $angle=acos($inner_pro/$dist1/$dist2);
    $sum_angle+=$angle*180/pi;
#    print "Angle is $sum_angle\n";
}
print $sum_angle/6,"\n";


