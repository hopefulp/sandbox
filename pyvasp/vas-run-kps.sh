#!/usr/bin/bash

##### setting parameter
x=3
N=4

##### default kpoints
kpoints=( 2 3 4 6 8 10 12 16 20 )
#echo ${kpoints[*]}

##### Usage
if [ $# -eq 0 ]; then
    exe=$(basename $0)
    echo "Usage: $exe POSCAR.name {kpoints list}"
    echo "Usage: $exe POSCAR.name ${kpoints[*]} "
    exit 1
fi

##### $1: POSCAR
poscar=$1
shift
fname=`cut -d. -f2 <<< $poscar`
echo $fname

kp=${@:-${kpoints[*]}}
#echo ${kp[*]}

for i in $kp; do
    echo $i
    echo "python $(which vas_make_ini.py) -s $poscar -j sp -d ${fname}k$i -i INCAR.sp -k $i $i 1 -x $x -N $N"
    done
