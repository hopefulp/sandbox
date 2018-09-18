#!/bin/bash

#dir1a=lda
#dir1b=relax2

dir1a=../vasp_relax/cpo27
dir1b=_CO2_done

dir2a=./Dos_6co2/DOScpo27
dir2b=_CO2_anal

for Me in `cat metal.dir`
  do
    dir=$dir1a$Me$dir1b
#    cp $dir/CONTCAR ./CONTCAR$Me
#    ./pos2xyz_park.pl CONTCAR$Me 2>/dev/null
    dir=$dir2a$Me$dir2b
#    cp $dir/ACF.dat ./ACF$Me.dat
#    echo $Me
#    ./xyz_dist.pl CONTCAR.xyz | awk '{print $4}'
    ./coulomb_mof72.pl CONTCAR$Me.xyz ACF$Me.dat $Me
  done


#mv $1$old $1$new
#mv $1$old.out $1$new.out

