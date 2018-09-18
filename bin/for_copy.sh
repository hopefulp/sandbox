#!/bin/bash

dir1=/qcfs/joonho/backup/PBE_D2_relax
base1=cpo27
base2=_cell_d2_r2
#base2=_CO2_r1_d2

dir2=~/vaspi/
#mkdir $dir2

for Me in `cat metal.dir`
    do
####	copy files into new directory

    	cp $dir1/$base1$Me$base2/CONTCAR $dir2/CONTCAR.MOF74.$Me
#	newdir=~/vaspi/$me
#	mkdir $newdir
#	cp $dir1/$base1$me$base2/CONTCAR $newdir

	#./xyz_dist.pl MOF_XYZ/$me.xyz
#	./xyz_angle.pl MOF_XYZ/$me.xyz
    done

