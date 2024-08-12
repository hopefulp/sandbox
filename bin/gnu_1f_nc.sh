#!/bin/bash

file=$1
shift

TITLE=$file

#yrange=$3

#term='postscript' 
#eps
#extens='png'
RANGE="[:] [:]"
if [ $# -eq 0 ]; then
    STYLE=" \"$file\" u 1:2  w l lt -1 lw 1 "
elif [ $# -ge 1 ]; then
    STYLE=" \"$file\" u 1:$1  w l lt 7 lw 1 "
fi
if [ $# -ge 2 ]; then
    STYLE="${STYLE}"", \"$file\" u 1:$2  w l lt 6 lw 1 "
fi
if [ $# -ge 3 ]; then
    STYLE="${STYLE}"", \"$file\" u 1:$3  w l lt 1 lw 1 "
fi
if [ $# -eq 4 ]; then
    STYLE="${STYLE}"", \"$file\" u 1:$4  w l lt -1 lw 1 "
fi

echo $STYLE

gnuplot -persist << EOF
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE  $STYLE
EOF

