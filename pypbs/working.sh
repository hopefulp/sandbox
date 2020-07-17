#!/bin/bash

if [ $# -eq 0 ]; then
    echo "$0 job( db fp1 tr )"
    exit 1
fi

job=$1

str=" -f OUTCAR -j tr -hl 4 -el 0.001 -fl 0.01 0.04 -nt 4000 -ntr 100 -dtype int -dl 1000 1001 -des gs -pf powNN -pn 5 -tef"
str_db=" -f ../OUTCAR -p amp-untrained-parameters.amp -im 1000 -ia 30 "
#tr="-db -i OUTCAR -qj NN5G2 -nc 10 -j tr -hl 4 -nt 4000 -ntr 100 -dtype int -dl 1000 1010 -m 3G -des gs -pf powNN -pn 5 -tef"
tr=" -f OUTCAR -j tr -hl 4 -el 0.001 -fl 0.01 0.04 -nt 4000 -ntr 10 -dtype int -dl 1000 1010 -des gs -pf powNN -pn 5 -tef -nc 10"

if [ $job == 'db' ]; then
    echo amp_run.py $str
    amp_run.py $str
elif [ $job == 'fp1' ]; then
    echo amp_anal.py $str_db
    amp_anal.py $str_db
elif [ $job == 'tr' ]; then
    echo amp_run.py $tr
    amp_run.py $tr
fi
