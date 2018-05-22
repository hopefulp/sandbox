#!/bin/tcsh
set pwdir = `pwd`

cd $pwdir
foreach sect (`ls -d */`)
	set fname = `basename $sect /`
	cd $fname

	if(-e INCAR) then
		if (`grep 'LWAVE = .TRUE.' INCAR|wc -l`) then
			sed -i -e 's|LWAVE = .TRUE.|LWAVE = .FALSE.|' INCAR
		else if (`grep 'LWAVE = .FALSE.' INCAR|wc -l`) then
			sed -i -e 's|LWAVE = .FALSE.|LWAVE = .TRUE.|' INCAR
		endif
		if (`grep 'LCHARG = .TRUE.' INCAR|wc -l`) then
			sed -i -e 's|LCHARG = .TRUE.|LCHARG = .FALSE.|' INCAR
		else if (`grep 'LCHARG = .FALSE.' INCAR|wc -l`) then
			sed -i -e 's|LCHARG = .FALSE.|LCHARG = .TRUE.|' INCAR
		endif
	endif
	cd $pwdir
end
