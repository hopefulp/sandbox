#!/usr/bin/perl
# written by Joonho Park
# read  DOSCAR 

$fin1=$ARGV[0];
$fin2=$ARGV[1];
$ref_sep=$ARGV[2];
$sys=$ARGV[3];


#system("./ldos_int.pl $fin1");
#system("./ldos_int.pl $fin2");

@energy1=();
@dos1=();
@energy2=();
@dos2=();

#open(NET, "netstat -i -n|") || die "can't fun netstat: $!";
open (IN, "./ldos_int.pl $fin1 $ref_sep |");
while (<IN>){
    chomp($_);
    @line=split(/\s+/,$_);
    if($line[0] eq ""){shift(@line);}
    push(@energy1,$line[0]);
    push(@dos1,$line[1]);
#    print $line[0],"\t",$line[1],"\n";
}
close IN; 

#open IN, '|-', "./ldos_int.pl $fin2"; 		# this prints in the screen
open (IN, "./ldos_int.pl $fin2 $ref_sep |");
while (<IN>){
    chomp;
    @line=split(/\s+/,$_);
    if($line[0] eq ""){shift(@line);}
    push(@energy2,$line[0]);
    push(@dos2,$line[1]);
#    print $line[0],"\t",$line[1],"\n";
}
close IN; 

if($#energy1>$#energy2) {$max_ndos=$#energy1;}
else		 	{$max_ndos=$#energy2;}

#for($i=0,$j=0,$k=0;$i<=$max_ndos;$i++){
$i=0,$j=0,$k=0;
if($sys=~/(M|m)e/){  $shift=0.5;}
else	          {  $shift=2;}  

while(){
	if(abs($energy1[$j]-$energy2[$k])<$shift){
	    printf "%10.3f  %15.5f  %10.3f  %15.5f %30.5f\n", $energy1[$j],$dos1[$j],$energy2[$k],$dos2[$k],
			$energy2[$k]*$dos2[$k]-$energy1[$j]*$dos1[$j];
	    $j++; $k++;
	}else{
	    if($energy1[$j]<$energy2[$k]){
		printf "%10.3f  %15.5f  %10.3f  %15.5f %30.5f\n", $energy1[$j],$dos1[$j], $energy1[$j] , 0.0,
			$energy2[$k]*0.0-$energy1[$j]*$dos1[$j];	
		$j++;
	    }else{
		printf "%10.3f  %15.5f  %10.3f  %15.5f %30.5f\n", $energy2[$j],0.0,$energy2[$k],$dos2[$k],
			$energy2[$k]*$dos2[$k]-$energy1[$j]*0.0;
		$k++;
	    }
	}
    	if($j>$#energy1 and $k>$#energy2) {
	    print "j= $j k= $k whereas energy1 $#energy1 and energy2 $#energy2\n";
	    last;
	}
#	    printf "%10.3f  %15.5f  %10.3f  %15.5f\n", $energy1[$i],$dos1[$i],$energy2[$i],$dos2[$i];
}

