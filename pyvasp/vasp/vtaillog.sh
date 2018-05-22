#!/bin/tcsh
set pwdir = `pwd`
set HOST = `hostname`

foreach sect (`ls -dv */`)
	set dirname = `basename $sect /`
	cd $dirname
	if (-e $dirname.log) then
		set type = 'log'
		set energy = `grep E0 *.log | tail -n 1 |awk '{print $5}'`
	else if (-e $dirname.out) then
		set type = 'out'
		set energy = `grep E0 *.out | tail -n 1 |awk '{print $5}'`
	else
		set type = ''
		set energy = ''
	endif
	echo $dirname '\t' $energy '\t' $type
	cd $pwdir
end
