#!/bin/bash

if [ $# != 3 ]; then
  echo "Input three dir's"
  exit
fi

infile1=$1
infile2=$2
infile3=$3
echo "Input file is $infile1 and $infile2"

dirname=`basename $PWD`
TITLE=$dirname
#TITLE="E-EF(eV)"

#term='postscript' 
#eps
term=png
extens=png
RANGE='[-15:] [:20]'
STYLE1=" u (\$1):(\$2)  w l lt 3 lw 0.5 "
STYLE2=" u (\$1):(\$2)  w l lt 1 lw 0.5 "
STYLE3=" u (\$1):(\$2)  w l lt -1 lw 0.5 "
gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$infile1"  $STYLE1, "$infile2" $STYLE2, "$infile3" $STYLE3
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

