#!/bin/bash

dir_pre=cpo27_1co2_
dir_lat=_dos

for Me in `cat metal.dir`
    do
	new_dir=$dir_pre$Me$dir_lat
    	~/bin/changeline.pl pbsfold-idft.csh $new_dir nodes=$1
    done
