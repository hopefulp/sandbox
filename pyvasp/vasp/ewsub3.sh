#!/bin/tcsh
set HOST = `hostname`

set onedir = `pwd`
foreach oned (`ls -dv */`)
	set onedname = `basename $oned /`
	cd $onedname

		set twodir = `pwd`
		foreach twod (`ls -dv */`)
			set twodname = `basename $twod /`
			cd $twodname

			set threedir = `pwd`
			foreach threed (`ls -dv con[0-1]*/`)
				set threedname = `basename $threed /`
				cd $threedname
				if(!(-e $threedname.out)&&!(-e $threedname.log)&&(-e INCAR)&&(-e POSCAR)&&(-e POTCAR)&&(-e KPOINTS)) then
					cpvasp.sh 1
					echo {$HOST} $threedname
					qsub vsp_{$HOST}.pbs
					sleep 1
				endif
				cd $threedir
			end
	
			cd $twodir
		end
	
	cd $onedir
end
