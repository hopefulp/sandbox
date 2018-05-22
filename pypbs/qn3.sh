#!/bin/tcsh

echo $HOST

if ($HOST == "rho") then
    set NCORE = 16
    set NNODES = 34
else if ($HOST == "psi") then
    set NCORE = 16
    set NNODES = 44
else if ($HOST == "kdft") then
    set NCORE = 12
    set NNODES = 24
else if ($HOST == "idft") then
    set NCORE = 8
    set NNODES = 22
else if ($HOST == "mu") then
    set NCORE = 32
    set NNODES = 10
else
    set NCORE = 1
    set NNODES = 50

    echo "hostname $HOST is not indentified"
    echo "default vaules for NCORE, NNODES $NCORE $NNODES"
endif

set IAM = `whoami`

echo "------------------------------"

set nfree = 0

foreach nds(01 02 03 04 05 06 07 08 09 )
	set ND = `qstat -nr | grep -o $HOST$nds/ | wc -l`
	set JOONHO = `qstat -nr | grep -B 1 $HOST$nds/ | grep $IAM | wc -l`
	set RS = `expr $NCORE - $ND`

	if ($RS == 0 && $JOONHO != 0) then
		echo "   [$nds] $ND / $NCORE ($IAM)"
	else if ($RS == 0 && $JOONHO == 0) then
		echo "   [$nds] $ND / $NCORE"
	else 
        echo "\033[33m * [$nds] $ND / $NCORE \033[32m($RS free procs) \033[37m"
        set nfree = `expr $nfree + 1 `
	endif
end

foreach nds(`seq 10 1 $NNODES`) 
        set ND = `qstat -nr | grep -o $HOST$nds/ | wc -l`
        set JOONHO = `qstat -nr | grep -B 1 $HOST$nds/ | grep $IAM | wc -l`
        set RS = `expr $NCORE - $ND`

        if ($RS == 0 && $JOONHO != 0) then
                echo "   [$nds] $ND / $NCORE ($IAM)"
        else if ($RS == 0 && $JOONHO == 0) then
                echo "   [$nds] $ND / $NCORE"
        else 
            echo "\033[33m * [$nds] $ND / $NCORE \033[32m($RS free procs) \033[37m"
            set nfree = `expr $nfree + 1 `
        endif
end
echo "------------------------------"
set NQ = `qstat -n | grep "Q   --" | wc -l`
echo "        $NQ queued jobs"
echo "        $nfree free nodes"
echo "------------------------------"
