#!/bin/tcsh
echo "------------------------------"

foreach nds(01 02 03 04 05 06 07 08 09)
	set ND = `qstat -nr | grep -o kdft$nds/ | wc -l`
	set JOONHO = `qstat -nr | grep -B 1 kdft$nds/ | grep joonho | wc -l`
	set RS = `expr 12 - $ND`

	if ($RS == 0 && $JOONHO != 0) then
		echo "   [$nds] $ND / 12 (joonho)"
	else if ($RS == 0 && $JOONHO == 0) then
		echo "   [$nds] $ND / 12"
	else echo "\033[33m * [$nds] $ND / 12 \033[32m($RS free procs) \033[37m"
	endif
end

foreach nds(10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25)
        set ND = `qstat -nr | grep -o kdft$nds/ | wc -l`
        set JOONHO = `qstat -nr | grep -B 1 kdft$nds/ | grep joonho | wc -l`
        set RS = `expr 12 - $ND`

        if ($RS == 0 && $JOONHO != 0) then
                echo "   [$nds] $ND / 12 (joonho)"
        else if ($RS == 0 && $JOONHO == 0) then
                echo "   [$nds] $ND / 12"
        else echo "\033[33m * [$nds] $ND / 12 \033[32m($RS free procs) \033[37m"
        endif
end

foreach nds(51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70)
        set ND = `qstat -nr | grep -o kdft$nds/ | wc -l`
        set JOONHO = `qstat -nr | grep -B 1 kdft$nds/ | grep joonho | wc -l`
        set RS = `expr 16 - $ND`

        if ($RS == 0 && $JOONHO != 0) then
                echo "   [$nds] $ND / 16 (joonho)"
        else if ($RS == 0 && $JOONHO == 0) then
                echo "   [$nds] $ND / 16"
        else echo "\033[33m * [$nds] $ND / 16 \033[32m($RS free procs) \033[37m"
        endif
end

echo "------------------------------"
set NQ = `qstat -n | grep "Q   --" | wc -l`
echo "        $NQ queued jobs"
echo "------------------------------"
