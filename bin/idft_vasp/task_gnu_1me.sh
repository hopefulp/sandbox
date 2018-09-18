#!/bin/bash

Me=$1

dir1a=cpo27_1co2_
dir2a=cpo27_1co2_ads_
dirb=_dos.extg

dir_des=$dir1a$Me$dirb
dir_ads=$dir2a$Me$dirb

cd $dir_des
    ../dosall.pl a=1 l=2 m=0 spin=

#./gnu_multi.sh $dir_des/Ldos_a1.dat $dir_des/Ldos_a25-51.dat $dir_ads/Ldos_a1.dat $dir_ads/Ldos_a25-51.dat
#./gnu_multi_arg.sh $dir_des/Ldos_a1.dat $dir_des/Ldos_a25-51.dat $dir_ads/Ldos_a1.dat $dir_ads/Ldos_a25-51.dat


