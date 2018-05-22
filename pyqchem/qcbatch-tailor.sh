#!/bin/bash



for x in `ls *-631G.mol`
  do
    f_pre1=$( basename $x | cut -d"-" -f1 )
    #f_pre2=$( basename $x | cut -d"-" -f2 )
    echo $f_pre1 #-$f_pre2
    echo "mv $x $f_pre1.mol"  # -$f_pre2.mol
  done

