#!/bin/bash

dir1a=cpo27_1co2_
dir1b=_new.d
dir2a=cpo27_1co2_ads_
dir2b=_dos.extg
#dir2b=_mag

#file=Ldos_a25-51.dat
file=Ldos_a1.dat
#file=Pdos_d_a1.dat
#file=Tdos.dat
for Me in `cat metal.dir`
  do
     dir_des=$dir1a$Me$dir1b
     dir_ads=$dir2a$Me$dir2b

     echo  $Me
#    ./get_vcharge.pl $dir_cell/ACF.dat | perl -ne 'chomp($_); print "$_\t";' ; ./get_vcharge.pl $dir_co2/ACF.dat
#      grep ZVAL $dir/POTCAR   

     cd $dir_des
#	../dosall.pl a=1 l=2
	../ldos_int_2fermi.pl $file 1.e-7 $Me des
     cd ../$dir_ads
	../ldos_int_2fermi.pl $file 1.e-7 $Me ads
     cd ..

  done
#	../dosall.pl a=1 l=2 s=1
#	../dosall.pl a=1 l=2 s=-1
#	../dosall.pl a=1 l=2 m=0 s=-1
#	../dosall.pl a=1 l=2 m=1 s=-1
#	../dosall.pl a=1 l=2 m=2 s=-1
#	../dosall.pl a=1 l=2 m=3 s=-1
#	../dosall.pl a=1 l=2 m=4 s=-1
#	../dosall.pl 
#	../dosall.pl a=1 l=0
#	../dosall.pl a=1 l=1
#	../dosall.pl a=1 l=2 
#        ../gnu1f.sh Tdos.dat
#        ../gnu1f.sh Pdos_d_a1.dat 

#    ./gnu_3f.sh $dir_cell/Pdos_da1-1.dat $dir_co2/Pdos_da1-1.dat $dir_co2/Pdos_pa32-32.dat
#     ./gnu_3f.sh $dir_des/Pdos_d_a1.dat $dir_des/Pdos_p_a25.dat $dir_ads/Pdos_d_a1.dat
#     ./gnu_3f.sh $dir_ads/Pdos_d_a1.dat $dir_des/Pdos_p_a25.dat $dir_ads/Pdos_p_a25.dat
#     ./gnu4f.sh $dir_des  $dir_ads Pdos_d_a1.dat Pdos_p_a25.dat $Me


