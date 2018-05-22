#!/bin/tcsh
foreach sect (`ls -dv 0?`)
	cat $sect/OUTCAR | grep 'energy  without entropy' | tail -n 1 | awk '{print $7}'
end

