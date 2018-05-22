#!/bin/tcsh
foreach sect (`ls *.msi`)
	msi2bgf.py $sect
end
