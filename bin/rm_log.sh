#!/bin/bash

if [ $# -lt 2 ]; then
    echo "it requires two index"
    echo "Usage:: $0 first_index last_index"
    exit
fi

ini=$1
fin=$2

i=$ini

while [ $i -le $fin ]
  do
    echo "rm $i.*"
    i=`expr $i + 1`
  done


