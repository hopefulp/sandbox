#!/bin/bash

for dir in `ls ac* -d`
    do
	echo $dir
	grep "Total energy" $dir/*out | awk '{print $10}' | perl -ne 'chomp($_); print "$_\t"; END {print "\n";}' | awk '{print $1 - $2 + $3 - $4}'

    done

