#!/bin/bash

prefix=MOF74-Fe-ccdc-1O2
molv=Acet4
mol=Ace

#for Me in `cat metal.txt`
#    do
#	fname=$prefix-$Me-$molv.pos
#	new_dir=$prefix-$Me-$molv
#	echo "file: $fname  dir: $new_dir"
#	./job_mkdir_sub_ff.sh $fname $new_dir $Me  $1 $mol
#    done

fname=$prefix.pos
new_dir=$prefix
./job_mkdir_sub_ff.sh $fname $new_dir Fe 1 1O2

