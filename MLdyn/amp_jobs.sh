#!/bin/bash
### in mlet
SB=/gpfs/home/joonho/sandboxg

source /home/joonho/bin/alias.sh

### job 0        1      2     3    4      5     for CASE
jobs=( "ampdir" "db" "wrapper1" "wrapper2" "te" "subdir" "sh" "chmod" )
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
outcar=$2
npart=$3    # sbatch in iron
idata=$4        # number of images per dir: fp 40 100
ndata=$5
dpre=NN     # G2p6max log200pn (des) NN (des) , part, fp (db), hl
dpre2=hl  # da  Nd300hl2020 Nd500hl2020 Nd500G34h202010         
dprehl=hl

#list=$(seq 11 1 20)                                # HL    , what is set?
#list=( 200 250 300 )
#list=( 7  9 12 14 16 18 20 )
#list=( 20 ) 
#list=(   7  9    12 14 )
#list=(  5  10 15  20 )
list=( 5 6 7 8 9 10  )
#list=(5)
#ndata=60
list2=$(seq $idata $ndata 900)                           # db
#list2=$(seq 60 $ndata 60)                           # db
#list2=($idata)
#echo $list2
#list=( 200 )
#list=$(seq 300 200 1000)
#list=( 200  800 1500 1800 2500 )      # input Ndata
#list=( 1000 1200 1500 1800 2000 )                   # input Ndata
#list2=( 20 30 40 50 )
listhl=( 7  9 10 12)
#listhl=(5 7)

pwd=`pwd`
### scan and make list
#list=$(ls NN* -d)
#echo list is $list
amp_dir=( "amp-neighborlists.ampdb" "amp-fingerprint-primes.ampdb" "amp-fingerprints.ampdb" )
#source /gpfs/home/joonho/bin/alias.sh
#echo $j
case $j in
    ### VASP make DIR for basis for the given descriptors: run training as serial
    ### run in twice nested subdir 
    "ampdir")
        dpre=NN     # G2p6max log200pn (des) NN (des) , part, fp (db), hl
        #dpre2=des
        for n in ${list[@]}; do
            dname=$dpre$n
            if [ ! -d $dname ]; then
                echo "amp_mkdir.py $dname -w amp -j des"
            fi
            echo "cd $dname"

            i=$idata
            j=$(expr $i + $ndata)
            ### NpowN
            #amp_wrapper.py -js qsub -j tr -qn NN9p$n -dl $i $j &
            if [ $(hostname) == 'login' ]; then
                echo "amp_wrapper.py -js qsub -j tr -qn $dname -dl $i $j  &"
            else
                ### sf NN$n or log10
                #echo "amp_wrapper.py -f $outcar -js sbatch -j tr -qn $dname -p $npart -sf log10 -dl $i $j -t db  &"
                echo "amp_wrapper.py -f $outcar -js sbatch -j tr -qn $dname -p $npart -sf NN$n -dl $i $j -t db  &"
                #echo "amp_wrapper.py -f $outcar -js sbatch -j tr -qn $dname -p $npart -np 20 -sf NN$n -dl $i  $j -t db  &"
            fi
            echo "cd .."
        done
        ;;
    ### make db for descriptor and for Ndata: 2 for-loop
    'db')
        dpre=NN     # G2p6max log200pn (des) NN (des) , part, fp (db), hl
        dpre2=db
        for n in ${list[@]}; do
            dname=$dpre$n
            ### run 1st data section with 'fp' to make ampdb's
            #echo "amp_mkdir.py $dname -w amp -j des"
            echo "cd $dname"
            for db in ${list2[@]}; do
                dname2=$dpre2$db
                echo "amp_mkdir.py $dname2 -w amp -j db"
                echo "cd $dname2"
                if [ $(hostname) == 'login' ]; then
                    echo "amp_wrapper.py -js qsub -j tr -qn pn${n}nd${db} -sf logpn${n}max200 -dl $db $(expr $db + 800) 3500 3600 -t db & "
                else
                    ### sf: NN5 or log10
                    echo "amp_wrapper.py -f $outcar -js sbatch -j tr -qn pn${n}nd${db} -p $npart -sf NN$n -dl $db $(expr $db + $ndata) -t db &"
                fi
                echo "cd .."

            done
            echo "cd .."
        done
        ;;
    ### Modify INCAR for DIR
    "wrapper1")
        for n in ${listhl[@]}; do
            hls=$n$n$n
            hlin="$n $n $n"
            dname=$dprehl$hls
            datain="0 800 1000"
            if [ ! -d $dname ]; then
                echo "amp_mkdir.py $dname -w amp -j tr "    # tr, db
            else
                echo "Do not overwrite amp-log.txt in the same directory; so exit"
                exit
            fi
            #echo "cp $dname2/amp-untrained-parameters.amp $dname2c"
            echo "cd $dname"
        #for dname in $(ls -d */); do
            #if [ ! -d $dname ]; then
            #    #echo "rm -r $dname"
            #    echo "amp_mkdir.py $dname -w amp -j des "       # for making db 'des' or 'tr'
            #fi
            #echo "cd $dname"
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
            echo "amp_wrapper.py -f $outcar -js sbatch -j trte -qn NN5hl$hls -p $npart -hl $hlin -sf NN5 -nt 1000-dl $datain &"
            echo "cd $pwd"
        done
        ;;
    #"wrapper_subdir")
    ### this is for training, if there is directory, stop qsub
    "wrapper2")                 # for HL
        dpre=NN     # G2p6max log200pn (des) NN (des) , part, fp (db), hl
        dpre2=db
        for n in ${list[@]}; do
            dname=$dpre$n
            echo "cd $dname"
        
            for n in ${listhl[@]}; do
                hls=$n$n$n
                hlin="$n $n $n"
                dname=$dprehl$hls
                if [ ! -d $dname ]; then
                    echo "amp_mkdir.py $dname -w amp -j tr "    # tr, db
                else
                    echo "Do not overwrite amp-log.txt in the same directory; so exit"
                    exit
                fi
                #echo "cp $dname2/amp-untrained-parameters.amp $dname2c"
                echo "cd $dname"
                ### Training
                ## mlet
                #echo "amp_wrapper.py -js qsub -qn Npn${n} -sf NN${n} -hl 8 8 -dl 1000 2000 3500 3800 & "
                #echo "amp_wrapper.py -js qsub -qn ${n} -sf ${n} -hl 8 8 -dl 1000 1500 3500 3800 & "
                ## iron
                #echo "amp_wrapper.py -f $outcar -js sbatch -j trte -qn Log5hl${n}${n} -p $npart -hl $n $n -sf log5 -dl 0 500 550 &"
                echo "amp_wrapper.py -f $outcar -js sbatch -j trte -qn NN5hl${n}${n} -p $npart -hl $hlin -sf NN5 -dl 0 500 550 &"
                ### Training-continue
                #echo "amp_wrapper.py -js qsub -qn NNpn${n}c -sf NN${n} -dl 1000 1300 3500 3600 & "
                ### TEST
                #echo "rm test*"
                #echo "amp_wrapper.py -js qsub -j te -qn NNpn${n}te -sf NN${n} -dl 1000 1300 3500 3700 & "
                ### database
                #echo "amp_wrapper.py -js qsub -qn NNpn${n}db -t db -sf NN${n} -hl 4 4 -dl 1100 2000 3500 4000 & "
                echo "cd .."
            done
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

