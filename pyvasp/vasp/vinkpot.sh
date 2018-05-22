#!/bin/tcsh
foreach sect (`ls -dv */`)
		cp INCAR KPOINTS POTCAR $sect
end
