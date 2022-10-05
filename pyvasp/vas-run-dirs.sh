#!/usr/bin/bash

##### setting parameter
### Pt
#x=3
#N=4

##### default kpoints
series=( 1Li 1LiP 2Li 2LiP 3Li 3LiP 4Li 4LiP )
#echo ${series[*]}
fpre=MXNBs22L1P
dseries=${@:-${series[@]}}
##### Usage
#if [ $# -eq 0 ]; then
#    dseries=$@
#else
#    dseries=$series[@]
#fi

for x in ${dseries[@]}; do
    #dirname=${fpre}${x}
    poscar=${fpre}${x}
    #echo $dirname
    #echo "python $(which vas_make_cont.py) -j opt -d ${dirname} "
    #python $(which vas_make_cont.py) -j opt -d ${dirname}
    #echo "python $(which vas_make_ini.py) -s $poscar -al -r"
    python $(which vas_make_ini.py) -s $poscar -al -r -rd
    done
