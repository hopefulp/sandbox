#!/bin/bash

jobkinds=( 'amp' 'qchem' )
subjobs=( 'node' 'qsub' )

if [ $# -eq 0 ]; then
    job='amp'
else
    job=$1
    shift
    if [ $# -ne 0 ]; then
        subjob=$1
        shift
    fi
fi

echo job=$job, subjob=$subjob

case $job in
    ### AMP jobs here
    'amp')
        para=log10
        param=${1:-$para}

        if [ $param == 'log10' ]; then
            str_param=" -des gs -pf log10 -pmm 0.05 200 -pn 10 -pmod del"
        elif [ $param == 'NpowN' ]; then
            str_param=" -des gs -pf powNN -pn 5"
        fi
        ### force train
        str_node="amp_run.py -f OUTCAR -j tr -hl 4 -el 0.001 -fl 0.1 0.04 -nt 4000 -ntr 100 -dtype int -dl 1000 1100 -nc 10"
        ### no force train
        str_node="amp_run.py -f OUTCAR -j tr -hl 4 -el 0.001 -fl -0.1  -nt 4000 -ntr 100 -dtype int -dl 1000 1100 1200 -nc 10 "
        str_node_suf=" -g"  # suppress plot
        ### NODE job here
        if [ -z "$subjob" -o "$subjob" == 'node' ]; then
            echo "In node:: "
            echo $str_node
            str_node="$str_node $str_param $str_node_suf" 
            echo "no forc training, no plot, w. test"
            echo $str_node

            echo GA_run.py -qj csh -nhls 5 -hl 5 -nn 12 \&
        fi
        ### Qsub job here
        if [ -z "$subjob" -o "$subjob" == 'qsub' ]; then
            echo -e "\n"
            echo "To get qsub command::"
            echo "sge_amp.py -db -i OUTCAR -qj g2offset -nc 10 -j tr -hl 4 -nt 4000 -ntr 100 -dtype int -dl 1000 1100 -m 3G -tef"
            echo "sge_amp.py -db -i OUTCAR -qj db -nc 10 -j tr -hl 4 -nt 4000 -ntr 1000 -dtype int -dl 0 1000 -m 9G -tef"
            echo "sge_amp.py -qj l200d12 -m 12G -hl 10 10 -el 0.001 -fl 0.1 -tef -nt 4000 -ntr 1000 -dtype div -dl 3 0 $str_param"

            echo -e "\nQsub command::"
            echo "qsub -N g2offset -pe numa 10 -l mem=3G -v fname=OUTCAR -v pyjob=tr -v hl=4 -v el=0.001 -v tef=force -v fl='0.1 0.04' -v nt=4000 -v ntr=100 -v dtype=int -v dlist='1000 1100' $SB/pypbs/sge_amp_module.csh"
            echo "qsub -N l200d12 -pe numa 12 -l mem=12G -v fname=OUTCAR -v pyjob=tr -v hl='10 10' -v el=0.001 -v tef=force -v des=gs -v pf=log10 -v pmod=del -v pmm='0.05 200.0' -v pn=10 -v fl='0.1' -v nt=4000 -v ntr=1000 -v dtype=div -v dlist='3 0' $SB/pypbs/sge_amp.csh"
            echo GA_run.py  -nhls 5 -hl 5 -nn 12
        fi
        ;;
    *)
        echo "input job in ( ${jobkinds[@]} )"
        ;;
esac        
