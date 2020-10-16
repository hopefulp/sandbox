#!/usr/bin/csh
#$ -cwd
#$ -N amp
#$ -V

## bash is not working
## PBS_VARIABLE is not working

if ( ! $?fname || ! $?pyjob || ! $?hl || ! $?fl ) then
    echo "variable fname, pyjob is not defined"
    echo "use:: qsub -N jobname -pe numa $nc -v fname=a.extxyz(OUTCAR) -v pyjob=tr $SB/pypbs/sge_amp.sh"
    exit(1)
endif    

set logfile = "$SGE_O_WORKDIR/sge$JOB_ID"
set nproc = $NSLOTS
set sandbox = "$HOME/sandboxg"

setenv PYTHONPATH $sandbox/pycommon:$sandbox/myplot:$sandbox/acpype:$sandbox/py_ai:$sandbox/chem_mod
set PYTHON = "$HOME/anaconda3/bin/python"
set EXE = "$sandbox/pyamp/amp_run.py"
#set EXE = "$sandbox/pyamp/amp_run_newest.py"

if ( ! $?nt ) then
    set nt = 100
    set ntr = 100
endif

set string = "-inf $fname -j $pyjob -hl $hl -el $el -fl $fl -nc $nproc"
if ( $?dtype ) then
    set string = "$string -nt $nt -ntr $ntr -dtype $dtype -dl $dlist"
endif
### in case: -des = gs, -pf powNN; if -pf log10, -pmm "pmin pmax" is added
if ( $?des ) then
    set string = "$string -des $des -pf $pf -pmod $pmod -pn $pn"
    if ( $pf == "log10" ) then
        set string = "$string -pmm $pmm"
    endif
endif
if ( $?tef ) then
    set string = "$string -tef"
endif
# depricate graphic in qsub, do not write comment line following a sentence
set string = "$string -g"

date >> $logfile
echo "Job options: $string" >> $logfile
echo "HOSTNAME    JOB_NAME   NSLOTS(nproc)" >> $logfile
echo "$HOSTNAME       $JOB_NAME        $NSLOTS " >> $logfile
$PYTHON $EXE $string
date >> $logfile
