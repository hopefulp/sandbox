#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage:: $0 Metal"
    echo "        Charge is imported from Metal.chg"
    echo "        Amol2mol.pl M-#.Amol charge M"
    echo "        # of line in M.chg and # of M-#.Amol files are supposed to be same"
    exit
fi

M=$1
chg_file=$M.chg
exe="/qcfs/joonho/bin/Amol2mol.pl"
i=1
f=12

for chg in $(cat $chg_file)
  do
    echo $chg
    echo "$exe $M-$i.Amol $chg $M"
    i=`expr $i + 1`
  done


