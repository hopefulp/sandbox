#!/bin/bash

file=$1

TITLE=$2

if [ $# -lt 2 ]; then
    TITLE=`basename $PWD`
fi


#term='postscript' 
#eps
#extens='png'
RANGE='[0:] [:]'
#STYLE1=" u ((\$1-$asymp)*$hkj)  w lp lt 1 lw 3 "
#STYLE1=" u (\$1):(\$8)  w l lt 1 lw 1 "
STYLE1=" u (\$8)  w l lt 1 lw 1 "
STYLE2=" u 1 w lp lt 29 lw 3 "
gnuplot -persist << EOF
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$file"  $STYLE1
EOF


