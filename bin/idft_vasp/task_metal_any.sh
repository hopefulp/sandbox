#!/bin/bash

dir1a=cpo27_1co2_
dir2a=cpo27_1co2_ads_
dirb=_dos

for Me in `cat metal.dir`
  do
     dir_des=$dir1a$Me$dirb
     dir_ads=$dir2a$Me$dirb

     ./get_fermi.pl $dir_des $dir_ads

#     cd $dir_des
#	../dosall1.pl l=2 a:1
#	../dosall1.pl l=1 a:25

 #    cd ../$dir_ads
#	../dosall1.pl l=2 a:1
#	../dosall1.pl l=1 a:25

  #  cd ..

  done
#     dir=$dir1$Me$dir2
#     cd $dir
#        ../gnu1f.sh Tdos.dat
#        ../gnu_2f.sh Pdos_d_a1.dat Pdos_p_a25.dat

#    ./gnu_2f.sh $dir_des/Tdos.dat $dir_ads/Tdos.dat
#     ./gnu_3f.sh $dir_ads/Pdos_d_a1.dat $dir_des/Pdos_p_a25.dat $dir_ads/Pdos_p_a25.dat
#     ./gnu4f.sh $dir_des  $dir_ads Pdos_d_a1.dat Pdos_p_a25.dat $Me
