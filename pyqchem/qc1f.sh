#!/bin/bash

outf=$1
rem=pbed_optsp

pre_f=`basename $1 | cut -d'.' -f1`


echo $input


cat $pre_f.mol rem.$rem > $pre_f.in
changeline_pbs_kdft_QC.pl pbs_qchem1.sh $pre_f

