#!/usr/bin/perl
# written by Joonho Park

$fin=$ARGV[0];
$Me=$ARGV[1];

open(IN,"<$fin");

$i=0;
while($line=<IN>){
   $energy=$line;
}
close(IN);

@field=split(/\s+/,$energy);
if($field[0] eq ""){
    shift(@field);
}

print "$field[2]  $Me\n";

