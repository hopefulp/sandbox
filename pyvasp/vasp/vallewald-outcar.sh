#!/bin/tcsh
set pwdir = `pwd`

cd $pwdir
foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	if ((-e $fname.out)) then
		set ewe = `grep 'Ewald energy' OUTCAR |tail -n 1|awk '{print $5}'`
		printf "%-15s %15s\n" $fname $ewe
	endif
	cd $pwdir
end
