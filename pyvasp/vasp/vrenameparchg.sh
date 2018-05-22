#!/bin/tcsh
foreach sect (`ls PARCHG.*`)
	set str2 = `echo $sect | cut -d'.' -f2`
	set str3 = `echo $sect | cut -d'.' -f3`
	mv $sect b${str2}k${str3}CHGCAR
end
