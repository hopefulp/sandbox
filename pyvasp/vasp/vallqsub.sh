#!/bin/tcsh
set pwdir = `pwd`
set HOST = `hostname`
set NND = 1

if ($#argv == 1) then
	set NND = $argv[1]
endif

foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	if(!(-e $fname.out)&&!(-e $fname.log)&&(-e INCAR)&&(-e POSCAR)&&(-e POTCAR)&&(-e KPOINTS)) then
		cpvasp.sh $NND
		echo '\033[32m'$fname'\033[0m'
		qsub vsp_${HOST}.pbs
		sleep 1
	endif
	cd $pwdir
end
