#!/bin/bash
### in mlet
SB=/gpfs/home/joonho/sandboxg

### job 0        1      2     3    4      5     for CASE
jobs=( "fp" "wrapper" "wrapper_subdir" "te" "subdir" "sh" "chmod" )
if [ $# -eq 0 ]; then
    echo "Usage:: $0 job=[ ${jobs[@]} ] arg2=[filename|dprefix] "
    echo "\$1=job  fp  : finger print calculation"
    echo "        wrapper: run amp_wrapper.py for train/test or each one"
    echo "        te  : wrapper-te only, which might be \$2 later"
    echo "        sh  : run bash to treat several directories for post-amp analysis"
    exit
fi
### $1: job in job list
j=$1            
### $2: filename, dirname, 
key=$2        # Energy check: test_energy.dat, Force check: test_fstat_acc.txt
#ndat=$3        # number of images per dir: fp 40 100
dpre=NN     # G2p6max log200pn (des) NN (des) , part, fp (db), hl
dpre2=Nd500hl88  #. Nd300hl2020 Nd500hl2020 Nd500G34h202010         

#list=$(seq 11 1 20)                                # HL    , what is set?
#list=( 200 250 300 )
#list=( 7  9 12 14 16 18 20 )
#list=( 20 ) 
#list=(   7  9    12 14 )
list=(  5     10       15 20 )
#list2=$(seq 0 800 3999)                           # db
#list=( 200 )
#list=$(seq 300 200 1000)
#list=( 200  800 1500 1800 2500 )      # input Ndata
#list=( 1000 1200 1500 1800 2000 )                   # input Ndata
#list2=( 20 30 40 50 )

pwd=`pwd`
### scan and make list
#list=$(ls NN* -d)
#echo list is $list
amp_dir=( "amp-neighborlists.ampdb" "amp-fingerprint-primes.ampdb" "amp-fingerprints.ampdb" )
source /gpfs/home/joonho/bin/alias.sh
#echo $j
case $j in
    ### VASP make DIR for CONTINOUS JOB
    "fp")
        for num in ${list[@]}; do
            dname=${dpre}$num
            echo "make_dir.py $dname -w amp -j db"
            echo "cd $dname"
            ### NpowN
            #amp_wrapper.py -js qsub -j tr -qn NN9p$num -dl $num  $(expr $num + 360) &
            echo "amp_wrapper.py -js qsub -j tr -qn $dpre$num -dl $num  $(expr $num + 400)  &"
            echo "cd .."
        done
        ;;
    ### make db for descriptor and for Ndata
    'db')
        for n in ${list[@]}; do
            dname=$dpre$n
            #echo "make_dir.py $dname -w amp -j des"
            echo "cd $dname"
            for db in ${list2[@]}; do
                dname2=$dpre2$db
                #echo "make_dir.py $dname2 -w amp -j des"
                #echo "cd $dname2"
                #echo "amp_wrapper.py -js qsub -j tr -qn pn${n}nd${db} -sf logpn${n}max200 -dl $db $(expr $db + 800) 3500 3600 -t db & "
                #echo "cd .."

            done
            echo "cd .."
        done
        ;;
    ### Modify INCAR for DIR
    "wrapper")
        for n in ${list[@]}; do
            ### HL directory
            #dname=$dpre$n$n
            dname=$dpre$n
        #for dname in $(ls -d */); do
            if [ ! -d $dname ]; then
                #echo "rm -r $dname"
                echo "make_dir.py $dname -w amp -j des "       # for making db 'des' or 'tr'
            fi
            echo "cd $dname"
            ### (1) max_param(k) in eta-logspace, log10(0.05..k)
            #echo "amp_wrapper.py -js qsub -qn $dname -k $n &"
            ### (2) scan Ndata
            #echo "amp_wrapper.py -js qsub -qn $dname -hl 10 10 -dl 1000 $(expr 1000 + $n) 3500 3700  & "
            #echo "amp_wrapper.py -js qsub -qn $dname -hl 10 10 -dl $n $(expr $n + 100) 3500 3600 & "
            ### (3) scan HL
            #echo "amp_wrapper.py -js qsub -qn $dname -hl $n $n -s &"
            #echo "amp_wrapper.py -js qsub -qn $dname -hl $n $n $key &"
            ### (4) pn G2
            #echo "amp_wrapper.py -js qsub -qn log200pn${n} -sf logpn${n}max200 -hl 10 10 -dl 1000 1500 3500 3700 & "
            ###     pn G2 & G3
            #echo "amp_wrapper.py -js qsub -qn NNpn${n} -sf NN${n} -hl 8 8 -dl 1000 1100 3500 3600 & "
            ### this is for test of amp_wrapper
            #echo "amp_wrapper.py -js qsub -qn N$n -c &"
            #echo "rm -r db*"
            #echo "cd .."
            echo "cd $pwd"
        done
        ;;
    "wrapper_subdir")
        #echo ${list[@]}
        for n in ${list[@]}; do
            dname=$dpre$n
            #dname=$n
            #if [ ! -d $dname ]; then
            #    echo "mkdir $dname"
            #fi
            echo "cd $dname"
            dname2=$dpre2
            #dname2c=${dpre2}cont
            #echo "rm -r $dname2c"
            if [ ! -d $dname2 ]; then
                echo "make_dir.py $dname2 -w amp -j tr "    # tr, db
            fi
            #echo "cp $dname2/amp-untrained-parameters.amp $dname2c"
            echo "cd $dname2"
            ### Training
            echo "amp_wrapper.py -js qsub -qn Npn${n} -sf NN${n} -hl 8 8 -dl 1000 2000 3500 3800 & "
            #echo "amp_wrapper.py -js qsub -qn ${n} -sf ${n} -hl 8 8 -dl 1000 1500 3500 3800 & "
            ### Training-continue
            #echo "amp_wrapper.py -js qsub -qn NNpn${n}c -sf NN${n} -dl 1000 1300 3500 3600 & "
            ### TEST
            #echo "rm test*"
            #echo "amp_wrapper.py -js qsub -j te -qn NNpn${n}te -sf NN${n} -dl 1000 1300 3500 3700 & "
            ### database
            #echo "amp_wrapper.py -js qsub -qn NNpn${n}db -t db -sf NN${n} -hl 4 4 -dl 1100 2000 3500 4000 & "
            echo "cd $pwd"
        done
        ;;
    "clean")
        for n in ${list[@]}; do
            dname=$dpre$n
            cd $dname
            #dname2=$dpre2
            #cd $dname2
            clean.py -w amp pbs -y
            cd $pwd
        done
        ;;
    "rm")
        for n in ${list[@]}; do
            dname=$dpre$n
            cd $dname
            echo "rm -r pot*"
            cd $pwd
        done
        ;;
    ### QSUB in SGE for DIR
    'te')
        for dir in `ls G2* -d` ; do
            echo "cd $dir"
            echo "amp_wrapper.py -js qsub -j te -qn $dir &"
            echo "cd .."
        done
        ;;
    ### COPY to multiple DIR
    'subdir')
        #source /gpfs/home/joonho/bin/alias.sh
        for n in ${list[@]}; do
            dname=$dpre$n
            cd $dname
            echo $dname
            wl $ampd
            #for ampd in ${amp_dir[@]}; do
            #    #mkdir $ampd
            #    #mkdir $ampd/loose
            #    for db in ${list2[@]}; do
            #        dname2=$dpre2$db
            #        for file in $(ls $dname2/$ampd/loose/*); do
            #            f=`basename $file`
            #            ### test directory merge with link
            #            #ln -s $dname2/$ampd/loose/$f $ampd/loose/$f
            #            mv $dname2/$ampd/loose/$f $ampd/loose/
            #        done
            #    done
            #done
            cd $pwd
        done
        ;;
    'scan2')
        for scan1 in ${list[@]}; do
            for scan2 in ${list2[@]}; do
                dname=$dpre$scan1$dpre2$scan2
                echo "make_dir.py $dname -w amp -j des "
                echo "cd $dname"
                echo "amp_wrapper.py -js qsub -qn $dname -hl 10 10 -sf logm${scan1}${dpre2}${scan2} -dl 1000 1300 3500 3700  & "
                echo "cd .."
            done
        done            
        ;;
    'ampd')
        for n in ${list[@]}; do
            dir=$dpre$n
            #head -3 $dir/$fname
            echo $dir
            cd $dir
            wl $ampd
            cd $pwd
        done
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
    'grep')
        for dname in $(ls -dv */ ); do
            echo $dname
            grep "optimization " $dname/Ndata300/amp-log.txt -B 1 | awk '/T/ {print $8, $10}'
        done
        ;;
    *)
        ### clean in many directoreis
        echo "cd $d"
        cd $d
        echo "clean.py -w vasp -j rm"
        echo "cd .."
        cd ..
        ;;
esac

