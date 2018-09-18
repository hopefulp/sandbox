#!/bin/bash

i=$1
final=$2
n=$3



while [ $i -le $final ]
do
    
    qdel $i.idft
    i=$(($i + 1))
done
    
