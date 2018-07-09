#!/bin/tcsh
set pwdir = `pwd`
set HOST = `hostname`

foreach sect (`ls -d */`)
	set fname = `basename $sect /`
	cd $fname
	if((-e POTCAR.Z)&&!(-e POTCAR)) then
		echo $fname
		zcat POTCAR.Z > POTCAR
	endif
	cd $pwdir
end
