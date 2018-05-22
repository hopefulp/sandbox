#!/bin/bash

outf=$1
rem=pbed_spaugth
pre_f=`basename $1 | cut -d'.' -f1`

input=$pre_f-sp

echo $input

qgeo_mo.pl $pre_f.out > $pre_f-augdz.xyz
xyz22mol.pl $pre_f-augdz.xyz

cat $pre_f-augdz.mol rem.$rem > $input.in
changeline_pbs_kdft_QC.pl pbs_qchem1.sh $input

