#!/bin/bash

rem=b3lyp_G631sD3_eda

remfile=${1:-$rem}

#pbs=pbs_qchem1_$HOST.sh

ippn=12
ppn=${2:-$ippn}

for molfile in `ls *.mol`
  do
    fmol=$( basename $molfile | cut -d'.' -f1 )
    echo "cat $molfile  rem.$remfile > $fmol.in"
    cat $molfile  rem.$remfile > $fmol.in
#    $bin/changeline_pbs_QChem.pl $pbs $fmol  $ppn
  done

