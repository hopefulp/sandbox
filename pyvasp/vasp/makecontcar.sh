#!/bin/tcsh
set pwdir = `pwd`
set currdir = `basename $pwdir`

mkdir contcar_${currdir}

foreach sect (`ls -l | grep ^d | awk '{print $9}'`)
	set fname = `basename $sect /`
    set fname2 = `basename $sect / |sed 's/\.//g'`
#	if (-e $fname/$fname.out || -e $fname/$fname.1st.out) then
	if (-e $fname/CONTCAR) then
		cp $fname/CONTCAR ./contcar_${currdir}/${fname2}CONTCAR
	endif
end
