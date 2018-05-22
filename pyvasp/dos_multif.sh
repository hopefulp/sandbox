#!/bin/bash


#dir="DOScpo27V_CO2_anal"

dir=$1
firsta=1
seconda=2
cd $dir

../dosall1.pl  a:$firsta  l=0
../dosall1.pl  a:$firsta  l=1
../dosall1.pl  a:$firsta  l=2
#./dosall1.pl  a:1 l=1


cd ..

#./dosall1.pl $dir a:1 l=2


