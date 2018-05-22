#!/bin/tcsh
set pwdir = `pwd`
set curdir = `basename $pwdir`

cd $pwdir
printf "%-15s %15s %15s %10s\n" 'File name' 'Energy-1st' 'Energy-2nd' 'Vol-diff'
echo '----------------------------------------------------------'
foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	vcomparevolume.py |grep CONTCAR
#			set e1st = `grep E0 $fname.1st.out |tail -n 1 |awk '{print $5}'`
	cd $pwdir
end
