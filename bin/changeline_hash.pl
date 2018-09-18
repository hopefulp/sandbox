#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

$fin=$ARGV[0];
$key=$ARGV[1];
$value=$ARGV[2];
@fname=split(/\./,$fin);

$fout="$fname[0]"."a.inp";

open(IN,"<$fin");
open(OUT,">$fout");

@f_line=();
### read three atoms's coordinates
$i=0;
while($line=<IN>){
#    $i++;
    #print $line;
    @field=split(/\s+/,$line);
    if($field[0] eq "") {shift(@field);}
    if($field[0] eq $key){
	$field[1] = $value;
	$line = "    $field[0]  $field[1]\n";
	push(@f_line,$line);
    }else{
	push(@f_line,$line);
    }
}

print OUT @f_line;
close(IN); 
close(OUT);

system("mv $fout $fin");

