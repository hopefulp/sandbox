#!/bin/tcsh
set pwdir = `pwd`
set fname = `basename $pwdir`
echo $fname
mv $fname.out $fname.1st.out
cp POSCAR POSCAR.1st
cp CONTCAR POSCAR
