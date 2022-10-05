#!/usr/bin/bash

##### setting parameter
### Pt
#x=3
#N=4

##### default kpoints
series=( a c d f g i j l )
#echo ${series[*]}
fpre=MXs22L2cT1
dseries=${@:-${series[@]}}
##### Usage
#if [ $# -eq 0 ]; then
#    dseries=$@
#else
#    dseries=$series[@]
#fi

for x in ${dseries[@]}; do
    dirname=${fpre}${x}
    #echo $dirname
    #echo "python $(which vas_make_cont.py) -j opt -d ${dirname} "
    python $(which vas_make_cont.py) -j opt -d ${dirname}
    done
