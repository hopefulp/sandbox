#!/bin/bash

file=$1
shift
TITLE=$file

col1=$1
col2=$2

title1='Training set'
title2='Validation set'
if [ $col1 -eq 6 ]; then
    ylabel='E_{RMSE} [eV/atom)]'
elif [ $col1 -eq 11 ]; then
    ylabel='F [eV/A)]'
fi
#yrange=$3

#term='postscript' 
#eps
#extens='png'
RANGE="[:] [:]"
if [ $# -eq 0 ]; then
    STYLE=" \"$file\" u 2:6  w l lt -1 lw 1 "
elif [ $# -eq 2 ]; then
    STYLE=" \"$file\" u 2:$1  w l lt 2 lw 1 title \"$title1\", \"$file\" u 2:$2  w l lt 4 lw 1 title \"$title2\" "
fi

echo $STYLE

gnuplot -persist << EOF
set xlabel 'Epoch'
set ylabel "$ylabel"
set title "$TITLE"
plot $RANGE  $STYLE
EOF

