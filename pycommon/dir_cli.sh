#!/usr/bin/bash

f=$1

for d in $(ls -d re0*); do
    #echo $d
    ### REMOVE directory
    #echo "rm -rf $d"
    ### copy to multiple destinations
    echo "cp $f $d"             # copy incar to all directories
    ### clean in many directoreis
    #echo "cd $d"
    #cd $d
    #dir_clean_p2.py -w vasp -j rm
    #echo "cd .."
    #cd ..
    ### make vasp 2nd directory
    #echo $d
    #pre=${d:0:2}
    #lat=${d:2}
    #d_new=re0_$lat
    #d_new=${d}
    #echo $d_new
    #mkdir $d_new
    #echo "vmake_2nd.py $d_new -d $d -j hybrid"
    done

