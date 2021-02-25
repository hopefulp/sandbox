#!/bin/bash
### in mlet
SB=/gpfs/home/joonho/sandboxg

### job 0        1      2     3    4      5     for CASE
job=( "vmake" "incar" "qsub" "cp" "rm" "chmod" "ampwrap" )
if [ $# -eq 0 ]; then
    echo "list job[0]vmake [1]incar [2]qsub [3]cp [4]rm [5]chmod"
    exit
fi

re='^[0-9]+$'
if [[ $1 =~ $re ]]; then
    j=${job[$1]}
else
    j=$1
fi
#echo $j
dname=hl
match=${2:-$dname}
f=$match
shift
shift
#for d in $(ls -d ${match}*); do
for d in $(ls -d */); do
#for d in $@; do
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
        "ampwrap")
            dirname=$match$d$d
            echo "make_dir.py $dirname -w amp -j tr"
            echo "cd $dirname"
            #echo "amp_wrapper.py -js qsub -qn HL$d$d -hl $d $d &"
            echo "amp_wrapper.py -js qsub -qn $d$d -hl $d $d &"
            echo "cd .."
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
        ### REMOVE dir as whole
        'rmr')
            echo "rm -r $d$f"
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
            ### clean in many directoreis
            echo "cd $d"
            cd $d
            dir_clean_p2.py -w vasp -j rm
            echo "cd .."
            cd ..
            ;;
    esac
    done

