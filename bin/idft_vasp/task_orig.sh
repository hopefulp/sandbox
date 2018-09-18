#!/bin/bash

dir_pre=../vasp_relax/cpo27
dir_later=_done

dira=_cell
dirb=_CO2



for Me in `cat metal.dir`
  do
    dir1=$dir_pre${Me}$dira$dir_later
    dir2=$dir_pre${Me}$dirb$dir_later
    cp $dir1/CHGCAR DOScpo27${Me}${dira}_anal
    cp $dir2/CHGCAR DOScpo27${Me}${dirb}_anal
#    y=`echo $x| sed "s/vap/vas/"`
#    echo "mv $x ${x}.vas"
#     mv $x ${x}.vas
#    cd ..

  done


#mv $1$old $1$new
#mv $1$old.out $1$new.out

