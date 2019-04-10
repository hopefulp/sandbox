#!/usr/bin/tcsh
#$ -cwd
#$ -N amp
#$ -V

## bash is not working
## PBS_VARIABLE is not working

if ( ! $?fname || ! $?py_job ) then
    echo "variable fname, job is not defined"
    echo "use:: qsub -v fname=a.extxyz -v job=tr sge_amp.sh"
    exit(1)
endif    

setenv PYTHONPATH $HOME/sandbox_gl/pycommon:$HOME/sandbox_gl/mymplot:$HOME/sandbox_gl/acpype:$HOME/sandbox_gl/py_ai
set PYTHON = "$HOME/anaconda3/bin/python"
set EXE = "$HOME/sandbox_gl/py_ai/amp_ene.py"

$PYTHON $EXE $fname $py_job -hl 4 4 4 -el 0.001 -g

