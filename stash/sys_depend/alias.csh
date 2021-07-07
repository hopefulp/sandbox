#!/bin/tcsh
set host = `hostname`
set kisti_id = 'x1274pjh'

alias idft "ssh -Y idft.kaist.ac.kr"
alias jdft "ssh -Y jdft.kaist.ac.kr"
alias kdft "ssh -Y kdft.kaist.ac.kr"
alias qch "ssh -Y qch.kaist.ac.kr"
alias cerius "ssh -Y cerius.kaist.ac.kr"
alias a4 "ssh -Y altix450.kaist.ac.kr"
alias a3 "ssh -Y altix3764.kaist.ac.kr"
alias hive "ssh -Y hive.wag.caltech.edu"
alias matrix "ssh -Y matrix.wag.caltech.edu"
alias psi "ssh  -Y psi.kaist.ac.kr"
alias rho "ssh  -Y rho.kaist.ac.kr"
alias tachyon "ssh $kisti_id@150.183.175.101"

alias runenv 'source /qcfs/joonho/scripts/.tcshrc.runenv'
alias ll "ls -lv"
alias lv "ls -v"
alias lt "ls -ltr"

set exe1  = /qcfs/joonho/Exe
set qtmp = /qcfs/joonho/tmp
set qcvin = /qcfs/joonho/binvasp
set qcvi = /qcfs/joonho/VaspINI
set qcvasp = /qcfs/common/vasp-pp
set qcvn = /qcfs/joonho/Silica/Npp
set qcv3 = /qcfs/joonho/V_work3
set qcw = /qcfs/joonho/Work
set qcww = /qcfs/joonho/singleM_CO2/Amine_O2
set qcwi = /qcfs/joonho/singleM_CO2/INIQChem

alias gow "cd /qcfs/joonho/Work"
alias gowi "cd $qcwi"
alias goww "cd $qcww"

alias gov "cd /qcfs/joonho/V_work"
alias gov1 "cd /qcfs/joonho/V_work1"
alias gov2 "cd /qcfs/joonho/V_work2"
alias gov3 "cd /qcfs/joonho/V_work3"
alias gov4 "cd /qcfs/joonho/V_work4"
alias gov5 "cd /qcfs/joonho/V_work5"
alias govn "cd $qcvn"
alias govi "cd /qcfs/joonho/VaspINI"

alias gosm "cd /qcfs/joonho/singleM_CO2"
alias gos "cd /qcfs/joonho/scripts"
alias goj "cd /qcfs/joonho"
alias goem "cd /qcfs/emergency/Joonho"

alias gobin "cd /qcfs/joonho/bin"
alias govin "cd /qcfs/joonho/binvasp" 
alias govasp "cd /qcfs/common/vasp-pp"
alias goe "cd $exe1"
alias psl "ps aux | grep out | wc -l"
alias gokhj "cd /qcfs/biduri/mof-74/hydrocarbon/from-exp"
alias lst "ls /qcfs/joonho/tmp"
alias lsbin "ls /qcfs/joonho/bin"
alias lsj "ls /qcfs/joonho"
alias lss "ls /qcfs/joonho/scripts"
alias source_a "source /qcfs/joonho/scripts/.tcshrc.alias"
alias vi_a "vi /qcfs/joonho/scripts/.tcshrc.alias"
alias vikenv "vi /qcfs/joonho/scripts/.tcshrc.kdft.env"
alias sourcekenv "source /qcfs/joonho/scripts/.tcshrc.kdft.env"
alias lsvasp "ls /qcfs/common/vasp-pp/"
alias ls "ls --color=tty"

# home dir
set hmj = /qcfs/joonho
set hmbin = /qcfs/joonho/bin
set hms = /qcfs/joonho/scripts
set hmback = /qcfs/joonho/backup

# QChem
alias gtot "grep 'total energy'" 
alias gTot "grep 'Total energy'"
alias gfin " grep 'Final energy' "
alias gfina " grep 'Final energy' *.out"
alias gret "grep Thank"
alias grett "grep Thank *.out"
alias greo "grep 'Optimization Cycle'"
alias gTOT "grep 'TOTAL ENERGY'"

alias aq "/qcfs/isty2e/scripts/allqsub"
alias qn "qstat -n"
alias qj "qstat | grep joonho; qstat | grep joonho | grep  ' R '  | wc -l "
alias qjj "qstat | grep 'joonho\|emergency'; qstat | grep 'joonho\|emergency' | grep R | wc -l "
alias qn3 "/qcfs/joonho/scripts/qn3.sh"
alias qu "qstat -u biduri"
alias qall "id qu; jd qu; kd qu"
#alias qd "/qcfs/biduri/scripts/qd.sh"
alias shq "/qcfs/biduri/scripts/shq.sh"
alias shq3 "/qcfs/biduri/scripts/shq3.sh"

if ( $host == "psi" ) then
    alias mo "/opt/applic/molden/molden"
else 
    alias mo "/qcfs/biduri/program/molden/molden"
endif

alias xx "/qcfs/biduri/scripts/vtstscripts/xdat2xyz.pl ; /qcfs/biduri/scripts/vtstscripts/pos2xyz.pl CONTCAR; /qcfs/biduri/scripts/vtstscripts/pos2xyz.pl POSCAR; cat POSCAR.xyz > out.xyz; cat CONTCAR.xyz >> out.xyz"
alias gp "gnuplot"
#alias gf "grep FORCES: OUTCAR; grep 'g(Force)' OUTCAR"
alias gf "grep 'g(F)' *.log; grep 'g(F)' *.out"
alias gm "grep 'number of electron' OUTCAR"
alias llog "ll */*.log"
alias llout "ll */*.out"
#alias vmd "/home/qclab/common/vmd/bin/vmd"

alias gref "grep F= "
alias greff "$qcvin/greff.sh"
alias dumax "du --max-depth=1 -h"

alias vesta "/qcfs/biduri/program/VESTA-x86_64-3.1.4/VESTA"
