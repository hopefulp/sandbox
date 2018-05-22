#!/bin/bash
### the number of arguments is one of 1, 2, 4 in the order of job_kind, file_w_list, file_frag1ix, file_frag2ix
#pre_dir=/qcfs/joonho/work/

job="cp"		
sys="mol"
frag1="FeX-1"
frag2="c2"
string="Total energy"

##########  0  1    2     3    4   5    6     7   8    9       10         11
job_list=( cp grep grepa grepf mv pchg clean rename rm rmr dos_fermi msi2pos.pl ls )

if [ $# -lt 1 ]; then
    echo "Usage::$0 job"
    exit
fi

Job=${1:-$job}
#Sys=${2:-$sys}
#Frag1=${3:-$frag1}
#Frag2=${4:-$frag2}
#String=${5:-$string}

for sys in `cat $Sys.dat`  	#`ls $frag1 -d1`; `cat $Sys.txt`
    do
	### for file
 #   	odir=$sys-$Frag1-$Frag2
#	ofile=$odir.out
	w_dir=$2

	if [ ! -d $w_dir -a ! -f $ofile ]; then
	    echo "There is not working directory: $w_dir! and no working file $ofile"
	    continue
	fi
	case $Job in
	    "${job_list[0]}")
		#echo "cp $w_dir/CONTCAR 1$sys-ads.vas"
		cp $w_dir/CONTCAR PBED2-6$sys.vas
		;;
	    "${job_list[1]}")
		echo "grep"
		grep F= $w_dir.out
		#grep F= $w_dir.out | awk '{print $3}'	
		#grep TOTEN $w_dir/OUTCAR
		;;
	    "${job_list[2]}")
		if [ -f $w_dir.out ]; then
		    grep F= $w_dir.out
		elif [ -f $w_dir.log ]; then
		    echo " There is $w_dir.log "
		    grep F= $w_dir.log
		fi
		;;
	    "${job_list[3]}")
		grep F= $w_dir.out  > ienergy.dat
		./anal_energy.pl ienergy.dat $sys $frag2
		;;
	    "${job_list[4]}")
		#echo "mv pec-$sys.vas  pec-$sys.outcar"
		for f in $(ls $w_dir/PARCHG*)
		    do
		    	echo "mv $f $f.vas"
		    done
		;;
	    "${job_list[5]}")
		cd $dir

		;;
	    "${job_list[6]}")
		tgt=$(echo $sys | sed -e "s/Fe-X/FeX/")
	    	if [ -e $tgt ]; then
	    	    echo "there is already the target file"
		    exit
	    	fi
	    	echo "mv $sys $tgt"
		;;
	    "${job_list[7]}")
		#echo "rm $w_dir/CHG"
		#echo "rm $w_dir/CHGCAR"
		echo "rm ./1$sys.vas"
		;;
	    "${job_list[8]}")
                echo "rm -r $w_dir"
                ;;
	    "${job_list[9]}")
		cd $w_dir
		    ../dosall.pl .
		    cd ..
		;;
	    "${job_list[10]}")
		echo "msi2pos.pl $odir C H"
		;;
            "${job_list[11]}")
		echo $w_dir
		ls $w_dir
		;;
	    *)
	   	echo "the keyword should be one of \" ${job_list[@]} \" "
		echo "present setting is no-job frag1ix=$frag1 $Sys-$frag2 " 
		exit
		;; 
	esac
#    	cd $dir
#    	for x in `ls -1 PARCHG*  `
#  	    do
#     	    	echo  $x
#      	    	mv $x $x.vas
#  	    done

#    	cd ..
#	grep F= $dir.out
    done
