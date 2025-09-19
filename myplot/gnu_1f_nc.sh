#!/bin/bash

file=$1
shift

TITLE=$file
legend1='dft-1/2'
legend2='dft'
#yrange=$3

#term='postscript' 
#eps
#extens='png'
RANGE='[:][:-15]'
if [ $# -eq 0 ]; then
    STYLE=" \"$file\" u 1:2  w l lt -1 lw 1 "
elif [ $# -ge 1 ]; then
    STYLE=" \"$file\" u 1:$1  w l lt 7 lw 1 title \"$legend1\" "
fi
if [ $# -ge 2 ]; then
    STYLE="${STYLE}"", \"$file\" u 1:$2 smooth cspline  w l lt 6 lw 1 title \"$legend2\" "
fi
if [ $# -ge 3 ]; then
    STYLE="${STYLE}"", \"$file\" u 1:$3  smooth cspline w l lt 1 lw 1 title \"$legend3\" "
fi
if [ $# -eq 4 ]; then
    STYLE="${STYLE}"", \"$file\" u 1:$4  w l lt -1 lw 1 "
fi

echo $RANGE $STYLE

gnuplot -persist << EOF
set xlabel 'Reaction coordinates'
set ylabel 'E [eV]'
set title "$TITLE"
plot [:][:-15]  $STYLE
EOF

