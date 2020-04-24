#!/bin/bash
hostname=`hostname`

function ln_source {
    for ampf in $@; do
        if [ -e $ampf ]; then
            echo "checked $ampf"
        elif [ -e ../$ampf ]; then
            ln -s ../$ampf $ampf
            echo "made $ampf"
        else
            echo "there is not amp data here and ../"
            exit 2
        fi
        done
    }

ampjob=tr

if [ $# -lt 1 ]; then
    echo "enter root directory name"
    exit 1
else
    dir=$1
    echo "mkdir $dir and cd $dir"
    if [ ! -e $dir ]; then
        mkdir $dir
    fi
    cd $dir
fi

if [ $ampjob == "tr" ]; then
    ampfiles=( 'amp-fingerprint-primes.ampdb' 'amp-fingerprints.ampdb' 'amp-neighborlists.ampdb' 'OUTCAR' )
elif [ $ampjob == "te" ]; then
    ampfiles=( 'amp-fingerprint-primes.ampdb' 'amp-fingerprints.ampdb' 'amp-neighborlists.ampdb' 'OUTCAR' 'amp.amp' )
fi

ln_source ${ampfiles[@]}

hl=( 6 7 8 9 )

for n in ${hl[@]}; do
    sub_dir=TR$n$n
    echo "mkdir $sub_dir and cd $sub_dir"
    if [ ! -e $sub_dir ]; then
        mkdir $sub_dir
    fi
    cd $sub_dir
    ln_source ${ampfiles[@]} 
    amp_run.py -f OUTCAR -j tr -nt 4500 -dt div -dl 3 0 -nc 16 -hl $n $n -el 0.001 -fl 0.00
    cd ..
    done

