#!/bin/bash

for file in `ls *out`
  do
    echo $file
    dir=`basename $file .out`
    if [ -d $dir ]; then
        echo "mv $file $dir"
    fi
  done
