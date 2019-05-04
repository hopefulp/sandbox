#!/bin/bash

cmd=$1


for f in $(pwd)/*; do
    fname=$(basename -- "$f")
    ext="${fname##*.}"
    fhead="${fname%.*}"
    if [ $ext == "out" ]; then
        echo "$cmd i=$fname"
    fi
done    
