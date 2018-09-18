#!/bin/bash

i=$1
final=$2
n=$3

host=`hostname`

while [ $i -le $final ]
do
    
    qdel $i.$host
    i=$(($i + 1))
done
    
