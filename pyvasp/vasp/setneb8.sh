#!/bin/tcsh
echo '' >> INCAR
echo '#neb' >> INCAR
echo 'IMAGES = 8' >> INCAR
echo 'SPRING = -5' >> INCAR
sed -i -e 's|#POTIM|POTIM|' INCAR
sed -i -e 's|#POTIM|POTIM|' INCAR
sed -i -e 's|ISIF = 3|ISIF = 2|' INCAR
sed -i -e 's|POTIM = 0.3|POTIM = 0.1|' INCAR
sed -i -e 's|POTIM = 0.2|POTIM = 0.1|' INCAR
sed -i -e 's|IBRION = [0-9]|IBRION = 3|' INCAR
#sed -i -e 's|EDIFF = 1E-5|EDIFF = 1E-4|' INCAR
sed -i -e 's|EDIFF = 5E-5|EDIFF = 1E-5|' INCAR
sed -i -e 's|EDIFFG = -0.025|EDIFFG = -0.1|' INCAR
sed -i -e 's|EDIFFG = 5E-4|EDIFFG = -0.1|' INCAR
sed -i -e 's|EDIFFG = 1E-3|EDIFFG = -0.1|' INCAR
sed -i -e 's|NPAR = [0-9]|NPAR = 2|' INCAR
sed -i -e 's|PREC = [a-zA-Z]*|PREC = normal|' INCAR

set curr_dir = `pwd`

cd $curr_dir
nebmake.pl POSCAR1 POSCAR2 8
cp $curr_dir/OUTCAR1 $curr_dir/00/OUTCAR
cp $curr_dir/POSCAR1 $curr_dir/00/CONTCAR
cp $curr_dir/OUTCAR2 $curr_dir/09/OUTCAR
cp $curr_dir/POSCAR2 $curr_dir/09/CONTCAR
