#!/bin/bash

ampjobs=( 'tr' 'te' 'log-plot' 'plot' )

function ln_source {
    for ampf in $@; do
        if [ -e $ampf ]; then
            echo "checked $ampf"
        elif [ -e ../$ampf ]; then
            ln -s ../$ampf $ampf
            echo "made $ampf"
        else
            echo "there is not amp data here and ../"
            exit 20
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

d_ampjob=tr
d_dname=Scan
d_subdir=TR

if [ $# -lt 1 ]; then
    echo "enter ampjob name"
    exit 1
else
    ampjob=${1:-$d_ampjob}
    shift
    if [ $ampjob == 'tr' ]; then
        if [ $# -lt 1 ]; then
            echo "enter root directory for scan for 'tr'"
            exit 2
        else
            dir=${1:-$d_dname}
            echo "mkdir $dir and cd $dir"
            if [ ! -e $dir ]; then
                mkdir $dir
            fi
            cd $dir
            shift
        fi
    fi    
fi
### this runs always
subdir=${1:-$d_subdir}

if [ $ampjob == "tr" ]; then
    ampfiles=( 'amp-fingerprint-primes.ampdb' 'amp-fingerprints.ampdb' 'amp-neighborlists.ampdb' 'OUTCAR' )
    ### for SCAN directory at first training
    echo "In scanning directory, call ln_source"
    #echo PWD `pwd | xargs basename`
    ln_source ${ampfiles[@]}
elif [ $ampjob == "te" ]; then
    ampfiles=( 'amp-fingerprint-primes.ampdb' 'amp-fingerprints.ampdb' 'amp-neighborlists.ampdb' 'OUTCAR' 'amp.amp' )
fi

hl=( 6 7 8 9 )

for n in ${hl[@]}; do
    ### check script
    sub_dir=${subdir}$n$n
    echo "check $sub_dir and cd $sub_dir"
    if [ ! -e $sub_dir ]; then
        mkdir $sub_dir
    fi
    echo "cd $sub_dir"
    cd $sub_dir
    case $ampjob in
        "tr")
            echo "call ln_source in subdirectory $sub_dir"
            ln_source ${ampfiles[@]} 
            #amp_run.py -f OUTCAR -j $ampjob -nt 4500 -dt div -dl 3 0 -nc 16 -hl $n $n -el 0.001 -fl 0.00
            amp_run.py -f OUTCAR -j $ampjob -nt 4500 -dt div -dl 3 0 -nc 16 -hl $n $n -el 0.001 -fl 0.1 -fc 0.04
            ;;
        "te")
            checkhost
            if [ ! -e 'test' ]; then
                mkdir test
            fi
            cd test
            ln_source ${ampfiles[@]}
            amp_run.py -f OUTCAR -j $ampjob -nt 4500 5000 -ntr 1500 -dt int -dl 0 500 -nc 4 -hl $n $n -el 0.001 -fl 0.00
            cd ..
            ;;
        "log-plot")
            checkhost
            amp-plotconvergence amp-log.txt
            cp convergence.pdf ../convergence$n$n.pdf
            ;;
        *)
            echo "There is not $ampjob in ampjob"
            ;;
        esac
    cd ..
    done

