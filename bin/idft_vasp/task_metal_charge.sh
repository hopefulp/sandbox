#!/bin/bash

#dir1=DOScpo27
#dir2a=_cell_anal
#dir2b=_CO2_anal

#dir1a=cpo27_1co2_
#dir1b=_new
#dir2a=cpo27_1co2_ads_
#dir2b=_mag

dira=cpo27
#dir1b=CO2_sp
#dirb=cell_sp
dirb=_CO2_r1_ch
for Me in `cat metal.dir`
  do
#    dir_cell=$dir1$Me$dir2a
#    dir_co2=$dir1$Me$dir2b
#     dir_des=$dir1a$Me$dir1b
#     dir_ads=$dir2a$Me$dir2b
     dir=$dira$Me$dirb
     echo  $Me
#    ./get_vcharge.pl $dir_cell/ACF.dat | perl -ne 'chomp($_); print "$_\t";' ; ./get_vcharge.pl $dir_co2/ACF.dat
    

#    cd $dir_cell
     cd $dir

	chgsum.pl AECCAR0 AECCAR2
	bader CHGCAR -ref CHGCAR_sum

#    cd ../$dir_co2
#     cd ../$dir_ads

#	chgsum.pl AECCAR0 AECCAR2
#	bader CHGCAR -ref CHGCAR_sum

    cd ..

  done



