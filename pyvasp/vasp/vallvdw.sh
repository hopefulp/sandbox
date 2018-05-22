#!/bin/tcsh
set pwdir = `pwd`
set curdir = `basename $pwdir`

cd $pwdir
foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	set evdw = `grep 'vdW energy (eV)' OUTCAR | tail -n 1 | awk '{print $5}'`
	printf "%-20s %10s\n" $fname $evdw
	cd $pwdir
end
