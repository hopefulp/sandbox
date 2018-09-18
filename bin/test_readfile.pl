#!/usr/bin/perl

my $finp=$ARGV[0];
print "INFILE: $finp\n";


open(IN,"<$finp");
while($line=<IN>){

    chomp($line);
    @line=split(/\s+/,$line);
    if($line[0] eq ""){ shift(@line);}
    print $line[0],"\n";
}
close(IN);

