#!/bin/bash

file=$1

TITLE=$2

yrange=$3

#term='postscript' 
#eps
#extens='png'
RANGE="[0:1] [:$yrange]"
#STYLE1=" u ((\$1-$asymp)*$hkj)  w lp lt 1 lw 3 "
STYLE1=" u (\$1):(\$2)  w l lt 1 lw 1 "
STYLE2=" u 1 w lp lt 29 lw 3 "
gnuplot -persist << EOF
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$file"  $STYLE1
EOF
#set term $term
#set out "$1.$extens"
#mkdir -p $HOME/public_html/badgeo
#mv $2-$3.$extens $HOME/public_html/badgeo
#mkdir -p $HOME/public_html/Spin


#mv $2-$3.$extens $HOME/public_html/Spin
#scp $2-$3.$extens joonho@cola.kaist.ac.kr:~/Spin

#set style line 1 lt -1 lw 1.5 pt 7 ps 0.8
#set style line 2 lt 1 lw 1.5 pt 5 ps 0.8
#set style line 3 lt 13 lw 1.5 pt 9 ps 0.8
# XX="1:"
# SET1='$2';
# SET2='$3';
#PLOT="\"$DAT\" u $XX(($SET1-asym)*fac) w lp ls 1 t \"fit surface($ver)\", \"$DAT\" u $XX(($SET2-asym)*fac) w lp ls 2 t \"ab initio\""

