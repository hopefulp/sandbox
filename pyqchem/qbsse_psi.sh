#!/bin/bash
# written by J Park
if [ $# -lt 4 ]; then
    echo "Usage:: $0 filename natom_A natom_B pbs_file ppn job_number"
    echo "	Make .in file and use \"pbs_qchem1_in_psi.sh\" "
    exit
fi

base=$1
na=$2
nb=$3
pbs_file=$4
tmp_ppn=16
tmp_job_number=5

ppn=${5:-$tmp_ppn}
job_number=${6:-$tmp_job_number}
#echo "ppn = $ppn"

fname=$base.in
dir=$base

mkdir $dir

if [ $? != "0" ]; then
    echo "There is $dir directory already"
    exit
fi

cp $pbs_file $dir
cp $fname $dir 

cd $dir
    /qcfs/joonho/bin/bsse_2mol_psi.pl $base $na $nb     # base name is without ".in"
    i=1
    if [ $job_number -eq 4 ]; then
	rm a$fname
    fi
#    exit
    for x in `ls *in`
	do
	    echo $x
    	    /qcfs/joonho/bin/changeline_pbs_bsse.pl $pbs_file $x $ppn  # select g1 for 12 process g2 for 16
	    i=`expr $i + 1`	# or i=$(($i + 1))
	done
    cd ..


