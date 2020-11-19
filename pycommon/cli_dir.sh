#!/bin/bash
### in mlet
SB=/gpfs/home/joonho/sandboxg

### job 0        1      2     3    4      5     for CASE
job=( "gitpush" "incar" "qsub" "cp" "rm" "chmod" )
if [ $# -eq 0 ]; then
    echo "list: ${job[@]}"
    exit
fi

j=$1

message='more'
arg1=${2:-$message}


case $j in
    ### VASP make DIR for CONTINOUS JOB
    "gitpush")
        echo "git add . -A"
        git add . -A
        echo "git commit -m \"$arg1\""
        git commit -m "$arg1"
        echo "git push"
        git push
        ;;
    ### Modify INCAR for DIR
    "incar")
        #echo "sed -i -e 's/Auto/.FALSE./'" $d/INCAR
        #echo "sed -i -e 's/ISTART = 0/ISTART = 1/' -e 's/ICHARG = 2/ICHARG = 0/' $d/INCAR"
        num=${d:6:3}
        echo "sed -i -e 's/ENCUT = 450/ENCUT = $num/' $d/INCAR"
        ;;
    ### QSUB in SGE for DIR
    'qsub')
        echo "qsub -N $d -pe numa 16 -v np=16 -v dir=$d $SB/pypbs/sge_vasp.csh"
        ;;
    ### COPY to multiple DIR
    'cp')
        echo "cp $f $d/INCAR"
        #echo cp $d/CONTCAR $d/POSCAR
        ;;
    ### REMOVE directory
    'rm')
        echo "rm $d/$f"
        ;;
    ### chmod for dir and files
    'chmod')
        if [ -d $d ]; then
            CHMOD=755
        else
            CHMOD=644
        fi
        echo chmod $CHMOD $d
        ;;
    *)
        echo "'$j' is not in the job list: (${job[@]})"
        exit
        ;;
esac

