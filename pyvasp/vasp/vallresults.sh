#!/bin/tcsh
#heejin
set pwdir = `pwd`
set curdir = `basename $pwdir`

cd $pwdir
printf "%-20s %15s %15s %10s %10s %5s\n" 'File name' 'Energy-1st' 'Energy-2nd' 'V-diff%' 'Max.Force' 'Step'
echo '--------------------------------------------------------------------------------'
foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	if((-e POSCAR)&&(-e INCAR)&&(-e POTCAR)&&(-e KPOINTS)) then
		set diff = `vcomparevolume.py | tail -n 1 | awk '{print $6}'`
#nothing
		if(!(-e ${fname}.log)&&!(-e ${fname}.out)&&!(-e ${fname}.1st.out)) then
			set e1st = ''
			set e2nd = ''
			set diff = ''
			set ff = ''
			printf "%-20s %15s %15s %10s %10s\n" $fname $e1st $e2nd $diff $ff
#1st running
		else if((-e ${fname}.log)&&!(-e ${fname}.out)&&!(-e ${fname}.1st.out)) then
			if (!(`grep 'TEBEG' INCAR|wc -l` )||(`grep '#TEBEG' INCAR|wc -l`)) then
				set ff = `gf | tail -n 1 | awk '{print $2}'`
				set e1st = `grep E0 $fname.log |tail -n 1 |awk '{print $5}'`
			else
				set ff = '-'
				set e1st = `vmdaverage.py |grep 'E  ='|awk '{print $3}'`
			endif
			set nstep = `grep E0 ${fname}.log |wc -l`
			printf "%-20s \033[0;33m%15s\033[0;0m                 %10s %10s %5s\n" $fname $e1st $diff $ff $nstep
#1st fin
		else if(!(-e ${fname}.log)&&(-e ${fname}.out)&&!(-e ${fname}.1st.out)) then
			if (!(`grep 'TEBEG' INCAR|wc -l` )||(`grep '#TEBEG' INCAR|wc -l`)) then
				set ff = `gf | tail -n 1 | awk '{print $2}'`
				set e1st = `grep E0 $fname.out |tail -n 1 |awk '{print $5}'`
			else
				set ff = '-'
				set e1st = `vmdaverage.py |grep 'E  ='|awk '{print $3}'`
			endif
			printf "%-20s %15s                 %10s %10s\n" $fname $e1st $diff $ff
#warning
		else if((-e ${fname}.log)&&(-e ${fname}.out)&&!(-e ${fname}.1st.out)) then
			printf "%-20s \033[0;31m  WARNING: .out file will be overwritten\033[0;0m\n" $fname
#1st fin, 2nd ready
		else if(!(-e ${fname}.log)&&!(-e ${fname}.out)&&(-e ${fname}.1st.out)) then
			set e1st = `grep E0 $fname.1st.out |tail -n 1 |awk '{print $5}'`
			set ff = `gf | tail -n 1 | awk '{print $2}'`
			printf "%-20s \033[0;32m%15s\033[0;0m                 %10s %10s\n" $fname $e1st $diff $ff
#1st fin, 2nd running
		else if((-e ${fname}.log)&&!(-e ${fname}.out)&&(-e ${fname}.1st.out)) then
			set e1st = `grep E0 $fname.1st.out |tail -n 1 |awk '{print $5}'`
			set e2nd = `grep E0 $fname.log |tail -n 1 |awk '{print $5}'`
			set ff = `gf | tail -n 1 | awk '{print $2}'`
			set nstep = `grep E0 ${fname}.log |wc -l`
			printf "%-20s %15s \033[0;33m%15s\033[0;0m %10s %10s %5s\n" $fname $e1st $e2nd $diff $ff $nstep
#1st fin, 2nd fin
		else if(!(-e ${fname}.log)&&(-e ${fname}.out)&&(-e ${fname}.1st.out)) then
			set e1st = `grep E0 $fname.1st.out |tail -n 1 |awk '{print $5}'`
			set e2nd = `grep E0 $fname.out |tail -n 1 |awk '{print $5}'`
			set ff = `gf | tail -n 1 | awk '{print $2}'`
			printf "%-20s %15s %15s %10s %10s %5s\n" $fname $e1st $e2nd $diff $ff
#warning
		else if((-e ${fname}.log)&&(-e ${fname}.out)&&(-e ${fname}.1st.out)) then
			printf "%-20s \033[0;31m  WARNING: .out file will be overwritten\033[0;0m\n" $fname
		endif
	endif
	if((-e POSCAR1)&&(-e POSCAR2)&&(-e OUTCAR1)&&(-e OUTCAR2)) then
		if((-e ${fname}.log)) then
			set nstep = `grep E0 $fname.log |wc -l`
			printf "\33[0;33m%s (step %s)\033[0;0m\n" $fname $nstep
			nebef.pl
		else if((-e ${fname}.out)) then
			set ntime = `grep 'Total CPU time used' 01/OUTCAR |awk '{print $6}'`
			set nstep = `grep E0 $fname.out |wc -l`
			printf "%s (finished, %s steps, %s sec)\n" $fname $nstep $ntime
			nebef.pl
		endif
	endif
	cd $pwdir
end
