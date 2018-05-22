#!/bin/tcsh
set pwdir = `pwd`
set HOST = `hostname`

foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	if(!(-e INCAR)) then
		echo $fname
		makevasp.py opt normal
	endif
	cd $pwdir
end
