#!/bin/bash
### in mlet
SB=/gpfs/home/joonho/sandboxg

### job 0        1      2     3    4  for CASE
job=( "vmake" "incar" "qsub" "cp" "rm")
if [ $# -eq 0 ]; then
    echo "list job[0]vmake [1]incar [2]qsub [3]cp [4]rm"
    exit
fi
j=${job[$1]}
#echo $j
match=$2
f=$3

for d in $(ls -d ${match}*); do
    case $j in
        ### VASP make DIR for CONTINOUS JOB
        "vmake")
            #echo $d
            pre=${d:0:2}
            lat=${d:2}
            d_new=re0$lat
            #echo $d_new
            echo "python /gpfs/home/joonho/sandboxg/pyvasp/vmake_cont.py $d $d_new -j hybrid"
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
        *)
            ### clean in many directoreis
            echo "cd $d"
            cd $d
            dir_clean_p2.py -w vasp -j rm
            echo "cd .."
            cd ..
            ;;
    esac
    done
