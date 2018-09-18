#!/usr/bin/perl
## written by Joonho Park
## fgaussian.pl sigma file or value
#

$PI=3.14159;
### Usage
if($#ARGV == -1){
    print "Usage: $0 sigma [filename(.dat)|points of point]\n";
    exit (1);
}

$sigma=shift(@ARGV);

if($#ARGV == 0 and $ARGV[0] =~ /\.freq/){
    $inf=$ARGV[0];
    @points=();
    open(IN,"<$inf");
    $i=0;
    while($line=<IN>){
	@line=split(/\s+/, $line);
	for($j=0;$j<=$#line;$j++){
	    if($line[$j] ne ""){
		push(@points,$line[$j]);
	    }
	}
    }
}else{
    @points=@ARGV;
}

close(IN);
print "array=", join(" ",@points), "\n";

$interval_in=$points[$#points]-$points[0];
#$offset=$interval_in/5;
$offset=3*$sigma;
$xini=$points[0]-$offset;
$xfin=$points[$#points]+$offset;
$interval=$xfin-$xini;
$ngrid=200;
$dx=$interval/$ngrid;

print "from $xini to $xfin with $ngrid grid and dx=$dx\n";

$i=0;
for($i=0;$i<$ngrid+1;$i++){
    $xi=$xini+$i*$dx;
    $intensity=0.;
    for($j=0;$j<=$#points;$j++){
	$intensity+=f_gaussian($xi,$points[$j]);
    	push(@intensity,$intensity);
    }
    print "$xi $intensity\n";
}

#print join("\n",@intensity), "\n";

sub f_gaussian{
    my ($x,$mu)=@_;
    my $f;

    $f=1./($sigma*sqrt(2.*$PI))*exp(-0.5*($x-$mu)*($x-$mu)/($sigma*$sigma));

    return $f;
}

