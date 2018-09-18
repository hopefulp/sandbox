#!/bin/bash

if [ $# != 4 -a $# != 5  ]; then
  echo "Input two dir and two files and Metal"
  exit
fi
dir1=$1
dir2=$2
filea=$3
fileb=$4
Me=$5

file1=$dir1/$filea
file2=$dir1/$fileb
file3=$dir2/$filea
file4=$dir2/$fileb
echo "Input file is $file1, $file2, $file3, and $file4"

#dirname=`basename $PWD`
#TITLE=$dirname
TITLE="$Me (Abs-Des)"

#term='postscript' 
#eps
term=png
extens=png
RANGE='[-13:8] [:20]'
STYLE1=" u (\$1):(\$2)  w l lt 2 lw 0.5 "
STYLE2=" u (\$1):(\$2)  w l lt 4 lw 0.5 "
STYLE3=" u (\$1):(\$2)  w l lt 3 lw 0.5 "
STYLE4=" u (\$1):(\$2)  w l lt 1 lw 0.5 "
STYLE5=" u (\$1):(\$2)  w l lt 5 lw 0.5 "
gnuplot -persist << EOF
#set term $term
#set out "$1$3.$extens"
set xlabel 'E(eV)'
set ylabel 'DOS'
set title "$TITLE"
plot $RANGE "$file1"  $STYLE1, "$file2" $STYLE2, "$file3" $STYLE3, "$file4" $STYLE4
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

