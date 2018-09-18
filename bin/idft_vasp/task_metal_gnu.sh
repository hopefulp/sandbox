#!/bin/bash

#dir1=DOScpo27
#dir2a=_cell_anal
#dir2b=_CO2_anal

dir1a=cpo27_1co2_
dir1b=_new.d
dir2a=cpo27_1co2_ads_
dir2b=_dos.extg

#dir1=../vasp_relax/cpo27
#dir2a=_cell_done
#dir2b=_CO2_done

for Me in `cat metal1.dir`
  do
#    dir_cell=$dir1$Me$dir2a
#    dir_co2=$dir1$Me$dir2b

     dir_des=$dir1a$Me$dir1b
     dir_ads=$dir2a$Me$dir2b

#     echo  $Me
#     dir=$dir1$Me$dir2a
#     head  $dir/POSCAR
#     head -2 $dir/POTCAR
#     dir=$dir1$Me$dir2b
#     head  $dir/POSCAR
#     head -2 $dir/POTCAR

     cd $dir_des
#	rm *dat
#        ../gnu1f.sh Tdos.dat
#        ../gnu1f.sh Pdos_d_a1.dat 

#	chgsum.pl AECCAR0 AECCAR2
#	bader CHGCAR -ref CHGCAR_sum

     cd ../$dir_ads

#	chgsum.pl AECCAR0 AECCAR2
#	bader CHGCAR -ref CHGCAR_sum

#	rm *dat
#	../dosall1.pl l=2 a:1
#	../dosall1.pl l=1 a:32
#        ../gnu1f.sh Tdos.dat
	
#     dir=$dir1$Me$dir2
#     cd $dir
#	../dosall1.pl l=2 a:1
#	../dosall1.pl l=1 a:32
#	../dosall1.pl l=1 a:25
#        ../gnu1f.sh Tdos.dat
#        ../gnu_2f.sh Pdos_d_a1.dat Pdos_p_a25.dat
      cd ..
    f1=Ldos_a1.dat
    f2=Ldos_a25-51.dat
#    f1=Tdos.dat

    ./gnu_multi.sh $dir_des/$f1 $dir_des/$f2 $dir_ads/$f1 $dir_ads/$f2
#    ./gnu_multi_arg.sh $dir_des/$f1 $dir_des/$f2 $dir_ads/$f1 $dir_ads/$f2
#    ./gnu_3f.sh $dir_cell/Pdos_da1-1.dat $dir_co2/Pdos_da1-1.dat $dir_co2/Pdos_pa32-32.dat
#     ./gnu_3f.sh $dir_des/Pdos_d_a1.dat $dir_des/Pdos_p_a25.dat $dir_ads/Pdos_d_a1.dat
#     ./gnu_3f.sh $dir_ads/Pdos_d_a1.dat $dir_des/Pdos_p_a25.dat $dir_ads/Pdos_p_a25.dat
#     ./gnu4f.sh $dir_des  $dir_ads Pdos_d_a1.dat Pdos_p_a25.dat $Me
  done


#mv $1$old $1$new
#mv $1$old.out $1$new.out

