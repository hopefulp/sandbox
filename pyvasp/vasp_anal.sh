#!/bin/bash

job_list=( bader bader2 convasp dos bchg end )
vbin=/qcfs/joonho/binvasp

if [ $# -lt 2 ]; then
    echo "Usage:: $0 job_name dir_name"
    echo "Job_name:: bader bader2 convasp dos bchg end"
    exit
else
    for job in ${job_list[@]} 
	do
	    #echo $job
	    if [ $1 == $job ]; then
		break
	    fi
	    if [ $job == "end" ]; then
		echo "Job list:: bader bader2 convasp dos bchg"
		exit
	    fi
	done
    job=$1
    echo "first argument $job is shifted as job name"
    shift
fi

if [ $1 == "*" ]; then
    $@=`find * -maxdepth 0 -type d`
fi

echo $@

for dir in  $@
    do
	
      case $job in
	"${job_list[0]}")
	    echo "cd $dir"
	    cd $dir
	    	chgsum.pl AECCAR0 AECCAR2
		bader CHGCAR -ref CHGCAR_sum
		cd ..
	    ;;
	"${job_list[1]}")
	    if [ -d $dir ]; then
	    $vbin/vasp_chgcar_del.pl $dir
	    echo "cd $dir"
	    cd $dir
		if [ -f "ACF2.dat" ]; then
		    echo "Already there is ACF2.dat !"
		elif [ -e "AECCAR2" ]; then
		    chgsum.pl AECCAR0 AECCAR2
	 	    chgsplit.sh CHGCAR
	 	    bader cf1 -ref CHGCAR_sum
	 	    mv ACF.dat ACF1.dat
	 	    mv AVF.dat AVF1.dat
	 	    mv BCF.dat BCF1.dat
	 	    bader cf2 -ref CHGCAR_sum
	 	    mv ACF.dat ACF2.dat
	 	    mv AVF.dat AVF2.dat
	 	    mv BCF.dat BCF2.dat
	 	    cd ..
		else
		    echo "Error:: There is not AECCAR2"
		fi
 	    fi
	    ;;
	"${job_list[2]}")
	    echo "cd $dir"
	    cd $dir
		convasp -chgint CHGCAR > x.chg &
		cd ..
	    ;;
	"${job_list[3]}")
	    cd $w_dir
	  	../dosall.pl .
		cd .. 
	    ;;
	"${job_list[4]}")
	    more $w_dir/ACF.dat | awk '{ print $5 }'
	    ;;	
	*)
	    echo "the keyword should be one of \" ${job_list[@]} \" "
                echo "present setting is no-job in $@ "
                exit
                ;;
      esac
    done
