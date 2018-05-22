#!/bin/tcsh
set pwdir = `pwd`
set curdir = `basename $pwdir`

cd $pwdir
foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	if (-e $fname.out) then
		mv $fname.out $fname
		rm $fname.pos $fname.k $fname.in $fname.pot
	else if (-e $fname.log) then
		mv $fname.log $fname
		rm $fname.pos $fname.k $fname.in $fname.pot
	endif
end

