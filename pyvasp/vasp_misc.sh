#!/bin/bash

job_list=( find findrm head )

vbin=/qcfs/joonho/binvasp
if [ $# -lt 2 ]; then
    echo "Usage:: $0 job_name dir_name filenames"
    exit
else
#    for job in ${job_list[@]} 
#	do
#	    #echo $job
#	    if [ $1 == $job ]; then
#		break
#	    fi
#	    if [ $job == "end" ]; then
#		exit
#	    fi
#	done
    Job=$1
    shift
    echo "first argument \"$Job\" is shifted as job name"
#    shift
fi

dir=$1
shift
fulldir=$PWD/$dir

echo "directry is saved in $dir and shifted"
case $Job in
    "${job_list[0]}")
	for file  in $*
	  do
	    find $fulldir -name $file
	  done
      ;;
    "${job_list[1]}")
	for file  in $*
	  do
	    if [ $file == "POSCAR" -o $file == "CONTCAR" ]; then
		echo " Are you sure to remove $file"
		exit
	    fi
	    find $fulldir -name $file -exec rm -rf {} \;
	  done
      ;;

  "${job_list[3]}")
      head $dir/ACF2.dat
      ;;
  "${job_list[4]}")
      echo "cd $dir"
      cd $dir
  	convasp -chgint CHGCAR > x.chg &
  	cd ..
      ;;
  "${job_list[5]}")
      cd $w_dir
    	../dosall.pl .
  	cd .. 
      ;;
  "${job_list[6]}")
      more $w_dir/ACF.dat | awk '{ print $5 }'
      ;;	
  *)
      echo "the keyword should be one of \" ${job_list[@]} \" "
          echo "present setting is no-job in $@ "
          exit
          ;;
esac
