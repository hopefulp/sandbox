#!/bin/bash

#dir=STRUCT


dir=/home/qclab/vasp/potpaw_GGA

dir1=Pot
frag=ggapaw


for x in `cat metal.list`
  do
#    cd $dir

#    y=`echo $x| sed "s/vap/vas/"`
#    echo "mv $x ${x}.vas"
#     mv $x ${x}.vas

#        cat $frag$x.pot ${frag}O.pot ${frag}C.pot ${frag}H.pot > $frag${x}OCH.pot
      
#    cd ..

# COPY potential
    if [ -e $dir/$x/POTCAR ]; then
      cp $dir/$x/POTCAR $dir1/${frag}_$x.pot
    elif [ -e $dir/$x/POTCAR.Z ]; then
      zcat $dir/$x/POTCAR.Z > $dir1/${frag}_$x.pot
    fi

    cd $dir1
        cat ${frag}_$x.pot ${frag}_O.pot ${frag}_C.pot ${frag}_H.pot > ${frag}_${x}OCH.pot
    cd ..

  done


#mv $1$old $1$new
#mv $1$old.out $1$new.out

