#!/bin/bash

bin=/qcfs/joonho/bin

for outfile in `ls *.out`
  do
    fmol=$( basename $outfile | cut -d'.' -f1 )
    #cat rem.eda_pbed $f_pre.mol > $f_pre.in
    echo $outfile
    grep -i jobtype $outfile | grep -i opt
    if [ $? == 0 ] ; then
        echo "opt"
        $bin/qgeo_mo.pl $outfile
    else
        echo "sp"
        $bin/qc_get_georem.pl m=$outfile
        $bin/xyz22mol.pl $fmol.mol
    fi
  done

