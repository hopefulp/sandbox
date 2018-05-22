#!/bin/tcsh
set pwdir = `pwd`

cd $pwdir
foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
		echo $fname
		mv $fname.log $fname.1st.log
		cp POSCAR POSCAR.1st
		cp CONTCAR POSCAR
#		sed -i -e 's|EDIFF = 1E-4|EDIFF = 1E-5|' INCAR
#		sed -i -e 's|EDIFF = 8E-4|EDIFF = 1E-5|' INCAR
#		sed -i -e 's|EDIFFG = -0.05|EDIFFG = -0.025|' INCAR
#		sed -i -e 's|PREC = normal|PREC = high|' INCAR
#		sed -i -e 's|#LMAXMIX = 6|LMAXMIX = 4|' INCAR
#		sed -i -e 's|LMAXMIX = 6|LMAXMIX = 4|' INCAR
#		sed -i -e 's|POTIM = 0.1|#POTIM = 0.1|' INCAR
#		sed -i -e 's|1 1 2|2 2 4|' KPOINTS
	#	sed -i -e 's|POTIM = 0.3|POTIM = 0.1|' INCAR
	#	sed -i -e 's|#POTIM = 0.1|POTIM = 0.1|' INCAR
	#	sed -i -e 's|#POTIM = 0.3|POTIM = 0.1|' INCAR
	cd $pwdir
end
