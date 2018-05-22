#!/bin/tcsh
set pwdir = `pwd`
set curdir = `basename $pwdir`

cd $pwdir
sed -i -e 's|ISIF = 3|ISIF = 2|' INCAR
cp CONTCAR templet
sed -i -e 's|   1.00000000000000|   Z.ZZ000000000000|' templet

#find INCAR -exec perl -pi -e "s|ISIF = 3|ISIF = 2|g" {} \;
#find templet -exec perl -pi -e "s|   1.00000000000000|   Z.ZZ000000000000|g" {} \;

foreach sect (0.80 0.85 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.00 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.10 1.15 1.20)
set wkdir = ${curdir}${sect}
mkdir $wkdir
find templet -exec perl -pi -e "s|   Z.ZZ|   $sect|g" {} \;
cp templet POSCAR
find templet -exec perl -pi -e "s|   $sect|   Z.ZZ|g" {} \;
cp INCAR POSCAR POTCAR KPOINTS ${pwdir}/${wkdir}
end

rm templet
rm POSCAR
