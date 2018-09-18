#!/bin/bash

for x in `cat metal.dir`
    do
#	./gnu2.sh DOScpo27${x}_cell_anal DOScpo27${x}_CO2_anal dos.dat
#	./gnu2f_pdos.sh DOScpo27${x}_cell_anal DOScpo27${x}_CO2_anal dosl2.dat
	./gnu2f_pdos.sh DOScpo27${x}_cell_anal DOScpo27${x}_CO2_anal Ldosa16.dat
    done
mv *png ../png
