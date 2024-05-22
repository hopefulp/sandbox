#!/bin/bash

TITLE="Potential"

#term='postscript' 
#eps
#extens='png'
RANGE='[2.5:7] [:0.25]'
STYLE1=" w lp lt 3 lw 1 "
STYLE2=" w lp lt 1 lw 1 "
STYLE3=" w lp lt 2 lw 1 "
STYLE4=" w lp lt 4 lw 1 "

gnuplot -persist << EOF
set xlabel 'r0'
set ylabel 'E(kcal/mol)'
set title "$TITLE"
D0O=0.1600E-00
r0O=3.4044

D0C=0.5589E-01
r0C=3.0946

DeO=0.1183
reO=3.7557
gammaO=10.1879

DeC=0.1270
reC=4.3625
gammaC=11.0671

plot $RANGE D0O*( (r0O/x)**12 - 2*(r0O/x)**6 ) $STYLE1,  DeO*( exp( -gammaO*(x/reO -1)) -2*exp( -(gammaO/2)*(x/reO -1))) $STYLE2,\
      D0C*( (r0C/x)**12 - 2*(r0C/x)**6 ) $STYLE3, DeC*( exp( -gammaC*(x/reC -1)) -2*exp( -(gammaC/2)*(x/reC -1))) $STYLE4




EOF
#, D0C*( (r0C/x)**12 - 2*(r0C/x)**6 ),  DeC*( exp( -gammaC*(x/reC -1)) -2*exp( -(gammaC/2)*(x/reC -1)))
#( (r0/x)^12 - 2*(r0/x)^6 ) 

