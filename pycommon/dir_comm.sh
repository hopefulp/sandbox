#!/usr/bin/bash

f=$1

for d in $(ls -d p*); do
    echo "cp $f $d"
    done

