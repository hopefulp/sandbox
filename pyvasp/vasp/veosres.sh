#!/bin/tcsh
set pwdir = `pwd`
set curdir = `basename $pwdir`
foreach sect (`ls -dv */`)
        set fname = `basename $sect /`
        set dis = `echo $fname | awk '{gsub(/'"${curdir}"'/,""); print}'`
        set en = `tail $sect/*out | grep E0 | awk '{print $5}'`
	echo $dis $en
#	echo $en
#       sed -n 2p ${fname}/energy
end
