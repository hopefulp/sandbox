#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage:: 1 dir_name and more than 2 f_name"
  exit
fi
infile1=$1/$2
infile2=$1/$3
infile3=$1/$4
infile4=$1/$5
infile5=$1/$6
infile6=$1/Fermi.dat
#echo "Input file is $infile2, $infile3, $infile3, $infile4, $infile5"

#dirname=`basename $PWD`
#TITLE=$dirname
TITLE="E(eV)"

#term='postscript' 
#eps
term=png
extens=png
RANGE='[:] [:5]'
STYLE1=" u (\$1):(\$2)  w l lt 2 lw 0.5 "
STYLE2=" u (\$1):(\$2)  w l lt 3 lw 0.5 "
STYLE3=" u (\$1):(\$2)  w l lt 1 lw 0.5 "
STYLE4=" u (\$1):(\$2)  w l lt 4 lw 0.5 "
STYLE5=" u (\$1):(\$2)  w l lt -1 lw 0.5 "
STYLE6=" u (\$1):(\$2)  w l lt -1 lw 1.0 "

if [ $# -eq 6 ]; then
gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$infile1"  $STYLE1, "$infile2" $STYLE2, "$infile3" $STYLE3, "$infile4" $STYLE4, "$infile5" $STYLE5,  "$infile6" $STYLE6
EOF

elif [ $# -eq 3 ]; then
gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$infile1"  $STYLE1, "$infile2" $STYLE2, "$infile6" $STYLE6
EOF

fi 
