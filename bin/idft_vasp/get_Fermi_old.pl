#!/usr/bin/perl
# written by Joonho Park

$fin=$ARGV[0];

$iatom_i=1;
$iatom_f=6;

open(IN,"<$fin");

$i=0;
$iatom=1;
$sum=0;
while($line=<IN>){
    next if($i<6);
    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){shift(@line);}
    next if($line[0] =~ /\D/);
#    print $line[0]."\n";
    if($iatom_i <= $line[0] and $line[0] <= $iatom_f){
	$sum+=$line[4];
    }
  
} continue {
    $i++;
}
close(IN);

print $sum/($iatom_f-$iatom_i+1),"\n";


