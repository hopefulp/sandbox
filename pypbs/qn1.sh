#!/bin/tcsh
echo "------------------------------"

foreach nds(1 2 3 4 5 6 7 8 9)
	set ND = `qstat -nr | grep -o node$nds/ | wc -l`
	set RS = `expr 8 - $ND`
	if ($RS != 0) then
		echo " * [ $nds] $ND / 8 ;\033[32m $RS free procs \033[37m"
		else echo "   [ $nds] $ND / 8" 
	endif
end

foreach nds(10 11 12 13 14 15 16 17 18 19 20)
	set ND = `qstat -nr | grep -o node$nds/ | wc -l`
	set RS = `expr 8 - $ND`
	if ($RS != 0) then
		echo " * [$nds] $ND / 8 ;\033[32m $RS free procs \033[37m"
		else echo "   [$nds] $ND / 8" 
	endif
end

echo "------------------------------"
set NQ = `qstat -n | grep "Q   --" | wc -l`
echo "        $NQ queued jobs"
echo "------------------------------"
