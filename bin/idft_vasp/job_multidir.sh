#!/bin/bash

echo $*

for dir in $@
  do
    cd $dir
	../dosall1.pl a:1 l=2
	../gnu1f.sh Tdos.dat
	../gnu1f.sh Pdos_d_a1.dat
    cd ..
  done

