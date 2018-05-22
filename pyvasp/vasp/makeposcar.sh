#!/bin/tcsh
set pwdir = `pwd`
set currdir = `basename $pwdir`

mkdir poscar_${currdir}

foreach sect (`ls -l | grep ^d | awk '{print $9}'`)
	set fname = `basename $sect /`
	if (-e $fname/POSCAR) then
		cp $fname/POSCAR ./poscar_${currdir}/${fname}POSCAR
	endif
#--new--
#	set fname = `basename $sect /`
#	set fname2 = `basename $sect / |sed 's/\.//g'`
#	cp $fname/POSCAR ./POSCAR_${currdir}/${fname2}POSCAR
#	cd POSCAR_${currdir}
#	convasp -direct < ${fname2}POSCAR > ${fname2}dPOSCAR
#	set line1 = `sed -n 1p ${fname2}POSCAR`
#	set line6 = `sed -n 6p ${fname2}POSCAR | awk '{print $1}'`
#	expr $line6 + 1 > /dev/null
#	if ($? == 0) then
#		sed -i 6" i\\$line1" ${fname2}dPOSCAR
#	endif
#	rm ${fname2}POSCAR
#	cd $pwdir
#	endif
end
