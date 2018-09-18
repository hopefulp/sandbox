#!/bin/bash

dir1=cpo27
dir2=_cell_d2

for Me in `cat metal.dir`
  do
#      echo $Me
      dir=$dir1$Me$dir2
      file=$dir.out
#      ./pos2xyz_park.pl $dir/CONTCAR
#      mv $dir/CONTCAR.xyz $Me.xyz
    
#     ./ldos_int_2fermi.pl  $dir1a$Me$dirb/Ldos_a25-51.dat 1.e-5 $Me
#      cp $dir/CONTCAR CONTCAR.$Me
#      ./get_lattice.pl cpo27${Me}cell_sp/CONTCAR
      #dir=/qcfs/joonho/MOF/${Me}_new.d
      #dir2=cpo27_1co2_${Me}_new.d
      #file2=$dir2.out

#      ./get_fermi.pl $dir_des $dir_ads
     grep F= $file | awk '{ print $3 }'  # ; grep F= $file2
#      echo "rm -r $dir $file"
#      rm -r $dir $file
#      cp -r $dir $dir2; cp $file $file2
     
#      ./vasp_ext_6meco2.pl $Me.xyz 
#     mv ${Me}_co2.mol ${Me}_co2.mol ../singleM_CO2/Me_co2   
#      grep ZVAL $dir/POTCAR   
#    mv ACF$Me.dat CONTCAR$Me* save
  done


#mv $1$old $1$new
#mv $1$old.out $1$new.out

