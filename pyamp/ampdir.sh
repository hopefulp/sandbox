#!/bin/bash
### in mlet
SB=/gpfs/home/joonho/sandboxg

### job 0        1      2     3    4      5     for CASE
jobs=( "fp" "wrapper" "te" "subdir" "sh" "chmod" )
if [ $# -eq 0 ]; then
    echo "Usage:: $0 job=[ ${jobs[@]} ] arg2=[filename|dir_prefix] "
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
dpre=eta     # G2p6max log200pn (des) NN (des) , part, fp (db), hl, hl1010nd 
dpre2=pn

#list=$(seq 11 1 20)                                 # HL    , what is set?
list=( 200 250 300 )
#list=( 5 40 80 100 200 250 300 )                   # fp log10
list2=$(seq 0 800 3999)                           # db
#list=( 200 )
#list=$(seq 300 200 1000)
#list=( 200  800 1500 1800 2500 )      # input Ndata
#list=( 1000 1200 1500 1800 2000 )                   # input Ndata
list2=( 20 30 40 50 )

amp_dir=( "amp-neighborlists.ampdb" "amp-fingerprint-primes.ampdb" "amp-fingerprints.ampdb" )
case $j in
    ### VASP make DIR for CONTINOUS JOB
    "fp")
        for num in ${list[@]}; do
            dname=${dir_pre}$num
            echo "make_dir.py $dname -w amp -j db"
            echo "cd $dname"
            ### NpowN
            #amp_wrapper.py -js qsub -j tr -qn NN9p$num -dl $num  $(expr $num + 360) &
            echo "amp_wrapper.py -js qsub -j tr -qn $dir_pre$num -dl $num  $(expr $num + 400)  &"
            echo "cd .."
        done
        ;;
    ### make db for descriptor and for Ndata
    'db')
        for n in ${list[@]}; do
            dname=$dir_pre$n
            #echo "make_dir.py $dname -w amp -j des"
            echo "cd $dname"
            for db in ${list2[@]}; do
                dname2=$dir_pre2$db
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
            #dname=$dir_pre$n$n
            dname=$dir_pre$n
            echo "make_dir.py $dname -w amp -j tr "       # for making db 'des' or 'tr'
        #for dname in $(ls -d */); do
            echo "cd $dname"

            ### (1) max_param(k) in eta-logspace, log10(0.05..k)
            #echo "amp_wrapper.py -js qsub -qn $dname -k $n &"
            ### (2) scan Ndata
            echo "amp_wrapper.py -js qsub -qn $dname -hl 10 10 -dl 1000 $(expr 1000 + $n) 3500 3700  & "
            #echo "amp_wrapper.py -js qsub -qn $dname -hl 10 10 -dl $n $(expr $n + 100) 3500 3600 & "
            ### (3) scan HL
            #echo "amp_wrapper.py -js qsub -qn $dname -hl $n $n -s &"
            #echo "amp_wrapper.py -js qsub -qn $dname -hl $n $n $key &"
            ### (4) pn
            #echo "amp_wrapper.py -js qsub -qn log200pn${n} -sf logpn${n}max200 -hl 10 10 -dl 1000 1500 3500 3700 & "
            ### this is for test of amp_wrapper
            #echo "amp_wrapper.py -js qsub -qn N$n -c &"
            #echo "rm -r db*"
            echo "cd .."
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
        for n in ${list[@]}; do
            dname=$dir_pre$n
            cd $dname
            for ampd in ${amp_dir[@]}; do
                #mkdir $ampd
                #mkdir $ampd/loose
                for db in ${list2[@]}; do
                    dname2=$dir_pre2$db
                    for file in $(ls $dname2/$ampd/loose/*); do
                        f=`basename $file`
                        ### test directory merge with link
                        #ln -s $dname2/$ampd/loose/$f $ampd/loose/$f
                        mv $dname2/$ampd/loose/$f $ampd/loose/
                    done
                done
            done
            cd ..
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
    'sh')
        for n in ${list[@]}; do
            dir=$dir_pre$n
            head -3 $dir/$fname
            cd $dir
            wl $ampd
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
    *)
        ### clean in many directoreis
        echo "cd $d"
        cd $d
        echo "clean.py -w vasp -j rm"
        echo "cd .."
        cd ..
        ;;
esac

