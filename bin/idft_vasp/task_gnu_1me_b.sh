#!/bin/bash

Me=$1
dir=$2
dir1a=cpo27_1co2_
dir2a=cpo27_1co2_ads_
dirb=_dos.extg

dir_des=$dir1a$Me$dirb
dir_ads=$dir2a$Me$dirb

#./gnu_multi.sh $dir_des/Pdos_s_a1.dat $dir_des/Ldos_a25-51.dat $dir_ads/Ldos_a1.dat $dir_ads/Ldos_a25-51.dat
./gnu_multi.sh $dir/Pdos_s_a1.dat $dir/Pdos_p_a1.dat  $dir/Pdos_d_a1.dat


