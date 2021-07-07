#!/bin/bash

ampd=( amp-fingerprints.ampdb amp-neighborlists.ampdb amp-fingerprint-primes.ampdb )

for dir1 in `ls -d */`; do
    for dir2 in ${ampd[@]}; do
        echo $dir1$dir2
        ls -al $dir1$dir2/loose | wc -l
        done
    done


    
