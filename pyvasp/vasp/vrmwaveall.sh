#!/bin/tcsh
set pwdir = `pwd`

cd $pwdir
foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	if (!(-e $fname.log)&&(-e $fname.1st.out)&&(-e $fname.out)) then
		if((-e WAVECAR)&&(-e CHG)) then
			rm WAVECAR CHG
			echo $fname '\t WAVECAR/CHG are removed'
		else
			echo $fname '\t WAVECAR/CHG are not found'
		endif
	endif
	cd $pwdir
end
