#!/bin/tcsh
set pwdir = `pwd`

foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	if(-e CONTCAR) then
		vclustercondel.py
	endif
	cd $pwdir
end
