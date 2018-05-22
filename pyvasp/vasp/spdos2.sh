#!/bin/tcsh
set pwd_dir = `pwd`
cd ..
set par_dir = `pwd`
cd $pwd_dir
#set curr_dir = `basename $par_dir`
set curr_dir = `basename $pwd_dir`

set dos_dir = ${curr_dir}-spdos11

mkdir $dos_dir
cp DOSCAR ./$dos_dir/
cp POSCAR ./$dos_dir/
cp CONTCAR ./$dos_dir/
cp OUTCAR ./$dos_dir/

cd $dos_dir
split_dos2

foreach sect (`ls DOS?`)
mv $sect $curr_dir-$sect
end

foreach sect (`ls DOS??`)
mv $sect $curr_dir-$sect
end

foreach sect (`ls DOS???`)
mv $sect $curr_dir-$sect
end

rm $curr_dir-DOSCAR
pos2xyz.pl POSCAR
head -n 6 ../DOSCAR > ./DOS
