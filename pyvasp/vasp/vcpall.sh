#!/bin/tcsh
set pwdir = `pwd`

foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	set newname = x${fname}
	vcp.sh $fname $newname
	cd $pwdir
end
