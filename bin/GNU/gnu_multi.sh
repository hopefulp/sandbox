#!/bin/bash

infile1=$1
infile2=$2
infile3=$3
infile4=$4
infile5=$5
infile6=$6

dirname=`basename $PWD`

TITLE="\"$dirname\" font \"Arial, 20\""
XLABLE="\"E(eV)\" font \" Helvetica, 15\""
YLABLE="\"DOS\" font \" Helvetica, 15\""
XTICS="font \"Times-Roman, 15\""
YTICS="font \"Times-Roman, 15\""
RANGE='[-15:][:20]'
KEY="font \", 8\""
#term='postscript' 
#eps
term=png
extens=png
STYLE1=" u (\$1):(\$2)  w l lt 3 lw 0.5 "
STYLE2=" u (\$1):(\$2)  w l lt 2 lw 1.0 "
STYLE3=" u (\$1):(\$2)  w l lt 1 lw 0.5 "
STYLE4=" u (\$1):(\$2)  w l lt -1  lw 0.5 "
STYLE5=" u (\$1):(\$2)  w l lt 4 lw 0.5 "
STYLE6=" u (\$1):(\$2)  w l lt 5 lw 0.5 "

if [ $# == 1 ]; then
gnuplot -persist << EOF
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$infile1"  $STYLE1
EOF

elif [ $# == 2 ]; then
gnuplot -persist << EOF
#set term $term
#set out "test.$extens"
set xlabel $XLABLE
set ylabel $YLABLE
set xtics $XTICS
set ytics $YTICS
set key font " ,8"
set title $TITLE
plot $RANGE "$infile1"  $STYLE3, "$infile2" $STYLE4
EOF

elif [ $# == 3 ]; then
gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set title $TITLE
set xlabel $XLABLE
set ylabel $YLABLE
set xtics $XTICS
set ytics $YTICS
set key font " ,8"
plot $RANGE "$infile1"  $STYLE1, "$infile2" $STYLE3, "$infile3" $STYLE4
EOF

elif [ $# == 4 ]; then
gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set title $TITLE
set xlabel $XLABLE
set ylabel $YLABLE
set xtics $XTICS
set ytics $YTICS
plot $RANGE "$infile1"  $STYLE1, "$infile2" $STYLE2, "$infile3" $STYLE3, "$infile4" $STYLE4
EOF

elif [ $# == 5 ]; then

gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set title $TITLE
set xlabel $XLABLE
set ylabel $YLABLE
set xtics $XTICS
set ytics $YTICS
plot $RANGE "$infile1"  $STYLE1, "$infile2" $STYLE2, "$infile3" $STYLE5, "$infile4" $STYLE3, "$infile5" $STYLE4
EOF

elif [ $# == 6 ]; then

gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$infile1"  $STYLE1, "$infile2" $STYLE2, "$infile3" $STYLE3, "$infile4" $STYLE4, "$infile5" $STYLE5, "$infile6" $STYLE6
EOF

fi

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

