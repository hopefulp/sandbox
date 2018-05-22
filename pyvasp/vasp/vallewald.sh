#!/bin/tcsh
set pwdir = `pwd`

cd $pwdir
if (-e ewaldout) then
	rm ewaldout
endif
touch ewaldout

foreach sect (`ls -dv */`)
	set fname = `basename $sect /`
	cd $fname
	cp CONTCAR temp
	sed -i 6d temp
	convasp -cart < temp > temp_cart
	convasp -names +2 +3.4 -2 +0.9 < temp_cart > temp_chg
#	convasp -names +2.8 +3.4 -2 +0.9 < temp_cart > temp_chg
	convasp -ewald < temp_chg |tail -n 1 >> $pwdir/ewaldout
	rm temp temp_cart temp_chg
	cd $pwdir
end
cat ewaldout
