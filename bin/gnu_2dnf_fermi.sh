#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Input 2 dir's and 1 filename"
  exit
fi

infile1=$1/$3
infile2=$2/$3
infile3=$1/Fermi.dat
infile4=$2/Fermi.dat
infile5=$1/$4
infile6=$2/$4
echo "Input file is $infile1 and $infile2"

#dirname=`basename $PWD`
#TITLE=$dirname
TITLE="E(eV)"

#term='postscript' 
#eps
term=png
extens=png
RANGE='[:] [:5]'
STYLE1=" u (\$1):(\$2)  w l lt 2 lw 0.5 "
STYLE2=" u (\$1):(\$2)  w l lt 1 lw 0.5 "
STYLE3=" u (\$1):(\$2)  w l lt 3 lw 0.5 "
STYLE4=" u (\$1):(\$2)  w l lt 4 lw 0.5 "

if [ $# -eq 3 ]; then
gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$infile1"  $STYLE1, "$infile2" $STYLE2, "$infile3" $STYLE1, "$infile4" $STYLE2
EOF

elif [ $# -eq 4 ]; then
gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$infile1"  $STYLE1, "$infile2" $STYLE2, "$infile3" $STYLE1, "$infile4" $STYLE2, "$infile5" $STYLE3, "$infile6" $STYLE4
EOF

fi
