#!/bin/bash

#dir1a=lda
#dir1a=gga
dir1=cpo27_1co2_
dir2=cpo27_1co2_ads_
#dir1b=relax2
#dir1b=co2_opt
#dir1a=../vasp_relax/cpo27
#dir1b=_CO2_done

#dir1=cpo27
#dir1b=cell_sp
#dir2b=CO2_sp

for Me in `cat metal.dir`
  do
#    dir=$dir1a$Me$dir1b
    dir=$dir1$Me
#    cp $dir/CONTCAR .
#    ./pos2xyz_park.pl CONTCAR 2>/dev/null
#    echo $Me
#    ./xyz_dist.pl CONTCAR.xyz | awk '{print $4}'

#    grep F= $dir.out > ienergy.dat
#    ./anal_energy.pl ienergy.dat $Me 
    grep F= $dir.out
#    dir=$dir1$Me$dir2b
    dir=$dir2$Me
    grep F= $dir.out

  done


#mv $1$old $1$new
#mv $1$old.out $1$new.out

