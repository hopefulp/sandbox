#!/bin/tcsh
set pwdir = `pwd`
set curdir = `basename $pwdir`

cd $pwdir
foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	echo compressing...  $fname
	gzip *
	cd $pwdir
end

