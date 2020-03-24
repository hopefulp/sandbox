#!/bin/bash

dir=`pwd`
i=0
def_k=0
k=${1:-$def_k}          # for case


bas2=GEN 
bas1="vdz(d)"
for fname in `ls $dir/*in`
    do
        f=`basename $fname`
        echo $f
        case $k in
            ### 1. change one line
            0)
                #echo "sed -i 's/$bas1/$bas2/'  $f"
                #sed -i 's/vdz(d)/GEN/' $f
                ### add one more line
                #echo "sed -i 's/$bas2/a PURECART 11'  $f"
                #sed -i '/GEN/a     PURECART     11' $f
                ### append file at the end of the other file
                cat bas.cc-pvdz >> $f
                ;;
            ### 2. queue submit all *.in file
            1)
                echo "Enter qjob_name nprocess nmem : for $f with "
                read qjob nproc nmem
                echo "qrun.sh $qjob $f $nproc $nmem: will you run? [Enter for yes]"
                read go
                if [ -z $go ]; then
                    qrun.sh $qjob $f $nproc $nmem
                    i=`expr $i + 1`
                else
                    exit $i
                fi
                ;;
        esac
    done
#echo $i job was loaded
