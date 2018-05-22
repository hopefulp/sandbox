#!/bin/tcsh
cp ../../$1/INCAR ./INCAR
cp ../../$1/KPOINTS ./KPOINTS
cp ../../$1/POTCAR ./POTCAR
cp ../../$1/CONTCAR ./POSCAR1
cp ../../$1/OUTCAR ./OUTCAR1
cp ../../$2/CONTCAR ./POSCAR2
cp ../../$2/OUTCAR ./OUTCAR2

#cp ../../geo/$1/INCAR ./INCAR
#cp ../../geo/$1/KPOINTS ./KPOINTS
#cp ../../geo/$1/POTCAR ./POTCAR
#cp ../../geo/$1/CONTCAR ./POSCAR1
#cp ../../geo/$1/OUTCAR ./OUTCAR1
#cp ../../geo/$2/CONTCAR ./POSCAR2
#cp ../../geo/$2/OUTCAR ./OUTCAR2
