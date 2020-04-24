#!/bin/bash

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

function checkhost {
    hostname=`hostname`
    if [ $hostname != "login" ]; then
        echo "run at login node"
        exit 10
    fi
    }

if [ $# -lt 1 ]; then
    echo "enter root directory name which has several sub_directories"
    exit 1
else
    dir=$1
    echo "mkdir $dir and cd $dir"
    if [ ! -e $dir ]; then
        mkdir $dir
    fi
    cd $dir
    shift
fi

ampjobs=( 'tr' 'te' 'log-plot' 'plot' )
d_ampjob=tr
if [ $# -lt 1 ]; then
    echo "there is no job"
    exit 2
else
    ampjob=${1:-$d_ampjob}
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
    case $ampjob in
        "tr")
            amp_run.py -f OUTCAR -j $ampjob -nt 4500 -dt div -dl 3 0 -nc 16 -hl $n $n -el 0.001 -fl 0.00
            ;;
        "te")
            checkhost
            amp_run.py -f OUTCAR -j $ampjob -nt 4500 5000 -ntr 1500 -dt int -dl 0 500 -nc 4 -hl $n $n -el 0.001 -fl 0.00
            ;;
        *)
            echo "Not ready except 'tr', 'te'"
            ;;
        esac
    cd ..
    done

