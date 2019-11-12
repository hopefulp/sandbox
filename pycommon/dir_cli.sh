#!/bin/bash
### in mlet
SB=/gpfs/home/joonho/sandboxg

match=$1

f=$2

for d in $(ls -d ${match}*); do
    #echo $d
    ### REMOVE directory
    #echo "rm -rf $d"
    ### copy to multiple destinations
    #echo "cp $f $d"             # copy incar to all directories, $SB/pyvasp/src/vdw_kernel.bindat
    #echo cp $d/CONTCAR $d/POSCAR

    ### Modify INCAR
    #echo "sed -i -e 's/Auto/.FALSE./'" $d/INCAR
    #echo "sed -i -e 's/ISTART = 0/ISTART = 1/' $d/INCAR"
    #echo "sed -i -e 's/ICHARG = 2/ICHARG = 0/' $d/INCAR"

    ### QSUB in SGE
    echo "qsub -N $d -pe numa 12 -v np=12 -v dir=$d $SB/pypbs/sge_vasp.csh"
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

