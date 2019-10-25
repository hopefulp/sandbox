#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
    echo "usage: $0 prefix datafile [2pt type = 1] [timestep] [grpfile|constrains] [rotsynumber=2] [linear=no]"
    exit(1)
endif

set prefix = $1
if (! -e ${prefix}.lammps) then
    echo "ERROR: Cannot find ${prefix}.lammps"
    exit(1)
endif
if (! -e ${prefix}.eng) then
    echo "ERROR: Cannot find ${prefix}.eng"
    exit(1)
endif
set mol = `basename $prefix`

set datafile = $2
if (! -e $2) then
    echo "ERROR: Cannot find data file $2"
    exit(1)
endif

set twoPt = 3
if ($#argv>2) then
    set twoPt = $3
endif
set memory = 2100
if (${twoPt} > 2) then
    set memory = 10500
    set mol = ${mol}.mol
endif

set tStep = 0.001
if ($#argv > 3) then
    set tStep = $4
endif

set grpString = ""
set fixString = ""
set rotSymString = "ANALYSIS_VAC_ROTN_SYMMETRY    2"
set linearString = "ANALYSIS_VAC_LINEAR_MOL	0"
if ($#argv > 4) then
    set fixString = "ANALYSIS_VAC_FIXED_DF        ${5}"
    if (-e $5) then
	set grpString = "IN_GROUPFILE                  ${5}"
        set fixString = "ANALYSIS_VAC_FIXED_DF         g"
	if($#argv < 5) then
	    set rotSymString = "ANALYSIS_VAC_ROTN_SYMMETRY    2"
	endif
	if($#argv < 6) then
	    set linearString = "ANALYSIS_VAC_LINEAR_MOL		g"
	endif
	set mol = ${mol}.grps
    endif
    if($#argv > 5) then
	set rotSymString = "ANALYSIS_VAC_ROTN_SYMMETRY    ${6}"
	if($#argv > 6) then
	    set linearString = "ANALYSIS_VAC_LINEAR_MOL ${7}"
	endif
    endif
endif

cat > ${mol}.in <<DATA;
IN_LMPDATA                    ${datafile}
IN_LMPTRJ                     ${prefix}
ANALYSIS_FRAME_INITIAL        1
ANALYSIS_FRAME_FINAL          0
ANALYSIS_FRAME_STEP           1
ANALYSIS_VAC_CORLENGTH        0.5
ANALYSIS_VAC_MEMORYMB         ${memory}
ANALYSIS_VAC_2PT              ${twoPt}
ANALYSIS_OUT                  ${mol}
ANALYSIS_LMP_TSTEP            ${tStep}
ANALYSIS_VAC_LINEAR_MOL	      0
${rotSymString}
${fixString}
${grpString}

DATA
