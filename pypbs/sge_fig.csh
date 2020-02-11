#!/usr/bin/csh
#$ -cwd
#$ -N amp
#$ -V

## bash is not working
## PBS_VARIABLE is not working

if ( ! $?fname || ! $?pyjob ) then
    echo "variable fname, pyjob is not defined"
    echo "use:: qsub -N jobname -pe numa $np -v fname=a.extxyz(OUTCAR) -v pyjob=tr -v np=$np $SB/pypbs/sge_amp.sh"
    exit(1)
endif    

set sandbox = "$HOME/sandboxg"

setenv PYTHONPATH $sandbox/pycommon:$sandbox/myplot:$sandbox/acpype:$sandbox/py_ai:$sandbox/chem_mod
set PYTHON = "$HOME/anaconda3/bin/python"
set EXE = "$sandbox/py_ai/amp_ene.py"

if ( ! $?scan ) then
    $PYTHON $EXE $fname $pyjob -tx -nc $np -hl 4 4 4 -el 0.0001 -g
else
    foreach i (`seq 2 2 10`)
        echo $i >> a.log
        ### repeat 5 times iteration
        foreach j  ( 1 2 3 4 5 )
            echo $i $j >> a.log
            ### validation check 1~4 parts among 5 parts
            if ( $pyjob == 'val' ) then
                foreach k  ( 0 1 2 3 )
                    echo $i $j $k >> a.log
                    $PYTHON $EXE $fname $pyjob -n 5 -hl $i $i $i -el 0.001 -i $k -g -nc $nc
                end
            else
                $PYTHON $EXE $fname $pyjob -n 5 -hl $i $i $i -el 0.001 -g
            endif
        end
    end
endif

