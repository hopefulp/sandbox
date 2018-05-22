#!/bin/bash

qcvin=/qcfs/joonho/binvasp
qcvi=/qcfs/joonho/VaspINI

w_dir=$1
tmp="full"
exe_vasp=${2:-$tmp}

#echo $exe_vasp

if [ $# -lt 1 ]; then
    echo "Usage:: $0 dir exe_vasp[full]"
    exit
fi

$qcvin/changeline_pbs_${HOST}_vasp.pl pbs-psi.csh $w_dir $exe_vasp
