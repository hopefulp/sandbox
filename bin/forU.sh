#!/bin/bash

U_I=1.5
U_diff=0.5
N_Upoint=19
N_TM=12
Mag="AFM"
job_list=( U_run OUTCAR gref )

if [ $# -lt 3 ]; then
    echo "Error in argument"
    echo "Usage:: $0 job dir1 dir2 \$4 \$5 \$6 n_TM AFM[FM]"
    echo "      with default variable of i_start= $U_I, interval= $U_diff, # of points= $N_Upoint"
    echo "	The job  should be one of \" ${job_list[@]} \" "
    exit 1
fi

qcvin=/qcfs/joonho/binvasp

job=$1
fname1=$2
fname2=$3

u_i=${4:-$U_I}
u_d=${5:-$U_diff}
n_upoint=${6:-$N_Upoint}
n_tm=${7:-$N_TM}
mag=${8:-$Mag}

u_f=`echo "$u_i + $u_d * ( $n_upoint - 1 ) " | bc`

echo " U = $u_i to $u_f" 

for i in $(seq 1 $n_upoint)
  do
    fnameout=${fname2}_u$u_i

    case $job in
      "${job_list[0]}")
	$qcvin/V_mkdir.tsh $fname1 $fnameout 1
	cd $fnameout
	  $qcvin/changeline_incar_U.pl U=$u_i
	  cd .. 
	$qcvin/changeline_pbs.pl pbs-psi-vasp.tsh $fnameout
	;;
      "${job_list[1]}")
	$qcvin/read_outcar.pl $fnameout 1 $n_tm $mag
	;;
      "${job_list[2]}")
	grep "F=" $fnameout.out
	;;
      *)
	echo "the keyword should be one of \" ${job_list[@]} \" "
	exit
	;;
    esac

    var=`echo "$u_i + $u_d" | bc `
    u_i=$var
  done
