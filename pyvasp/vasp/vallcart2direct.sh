#!/bin/tcsh
foreach sect (`ls *`)
	convasp -direct < $sect > dir$sect
end
