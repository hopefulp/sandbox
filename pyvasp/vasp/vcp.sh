#!/bin/tcsh
if ($#argv == 0) then
        echo "Usage : [source] [target]"
        exit()
else if ($#argv == 1) then
        echo "Usage : [source] [target]"
        exit()
else if ($#argv == 2) then
        set src = $argv[1]
        set tgt = $argv[2]
endif

if (!(-d $tgt)) then
	mkdir $tgt
	if (-e $src/INCAR.gz) cp $src/INCAR.gz $tgt/INCAR.gz
	if (-e $src/POSCAR.gz) cp $src/POSCAR.gz $tgt/POSCAR.gz
	if (-e $src/POTCAR.gz) cp $src/POTCAR.gz $tgt/POTCAR.gz
	if (-e $src/KPOINTS.gz) cp $src/KPOINTS.gz $tgt/KPOINTS.gz
	if (-e $src/CONTCAR.gz) cp $src/CONTCAR.gz $tgt/CONTCAR.gz
	if (-e $src/1stPOSCAR) cp $src/1stPOSCAR $tgt/1stPOSCAR
	if (-e $src/1POSCAR) cp $src/1POSCAR $tgt/1POSCAR
	cp $src/INCAR $tgt/INCAR
	cp $src/POSCAR $tgt/POSCAR
	cp $src/POTCAR $tgt/POTCAR
	cp $src/KPOINTS $tgt/KPOINTS
	cp $src/CONTCAR $tgt/CONTCAR
	cp $src/CONTCAR $tgt/POSCAR
else
	echo "$tgt already exist"
	exit()
endif
