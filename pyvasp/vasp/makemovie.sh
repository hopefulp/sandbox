#!/bin/tcsh
set curr_dir = `pwd`
rm $curr_dir/*.xyz

cp 00/POSCAR 00/CONTCAR
cp 09/POSCAR 09/CONTCAR

touch $curr_dir/movie.xyz
touch $curr_dir/poscarmovie.xyz
touch $curr_dir/contcarmovie.xyz

foreach sect (`ls -dv 0?`)
	cd $curr_dir/$sect
	rm *.xyz
	xx
	cat movie.xyz >> $curr_dir/movie.xyz
	cat POSCAR.xyz >> $curr_dir/poscarmovie.xyz
	cat CONTCAR.xyz >> $curr_dir/contcarmovie.xyz
end
