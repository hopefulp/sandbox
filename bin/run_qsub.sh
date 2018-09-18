#!/bin/bash

for dir in `cat odd_job.dir`
    do
	echo $dir
	~/bin/changeline.pl pbsfold-idft.csh ${dir} nodes=1
    done


