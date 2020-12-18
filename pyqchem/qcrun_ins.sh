#!/bin/bash

inode=1

i=0

for infile in `ls *.in`
  do
    fin=$( basename $infile | cut -d'.' -f1 )
    ### parallel
    #echo "qsub_server.py qchem -qj job$i -i ${fin}.in -n 16 -m 3 -no skylake@node$i"
    ### serial
    echo "qchem $fin.in $fin.out &"
    #qsub_server.py qchem -qj job${node_number} -i ${fin}.in -n 6 -m 3 -no skylake@node${$i}
#    $bin/changeline_pbs_QChem.pl $pbs $fmol  $ppn
  done

