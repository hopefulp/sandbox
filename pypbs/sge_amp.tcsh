#!/usr/bin/tcsh
#$ -cwd
#$ -N amp
#$ -pe numa 4
#$ -V

## bash is not working
## PBS_VARIABLE is not working

if ( ! $?fname || ! $?py_job ) then
    echo "variable fname, job is not defined"
    echo "use:: qsub -v fname=a.extxyz -v job=tr sge_amp.sh"
    exit(1)
endif    

/gpfs/home/joonho/sandbox_gl/py_ai/sge_amp_ene.py $fname $py_job -hl 4 4 4 -el 0.001 +g

