#!/bin/bash
#usage ./     filename 12 3 

base=$1
na=$2
nb=$3
pbs_file=$4
#group=$5
ppn=$5

if [ $ppn = "" ]; then
    ppn=12
fi

fname=$base.in
dir=$base

echo "ppn = $ppn"

# exit

mkdir $dir

if [ $? != "0" ]; then
    echo "There is $dir directory already"
    exit
fi

cp $pbs_file $dir
cp $fname $dir 
cd $dir
    ~/bin/bsse_2mol_paral.pl $base $na $nb     # base name is without ".in"
    i=1
    for x in `ls *in`
	do
	    echo $x
	    #if [ $i -eq 1 ]; then
	#	cat ../rem.rimp2_freq_add >> $x
	    #fi
    	    /qcfs/joonho/bin/changeline_pbs_bsse.pl $pbs_file $x $ppn  # select g1 for 12 process g2 for 16
	    #exit
	    i=`expr $i + 1`	# or i=$(($i + 1))
	done
    cd ..


