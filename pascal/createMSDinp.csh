#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
    echo "usage: $0 prefix datafile [timestep] [grpfile]"
    exit(1)
endif

set prefix = $1
if (! -e ${prefix}.lammps) then
    echo "ERROR: Cannot find ${prefix}.lammps"
    exit(1)
endif
set mol = `basename $prefix`

set datafile = $2
if (! -e $2) then
    echo "ERROR: Cannot find data file $2"
    exit(1)
endif

set tStep = 0.001
if ($#argv > 2) then
    set tStep = $3
endif

set grpString = ""
if ($#argv > 3) then
    set grpString = "IN_GROUPFILE                  ${4}"
endif

cat > ${mol}.in <<DATA;
IN_LMPDATA                    ${datafile}
IN_LMPTRJ                     ${prefix}
ANALYSIS_FRAME_INITIAL        1
ANALYSIS_FRAME_FINAL          0
ANALYSIS_FRAME_STEP           1
ANALYSIS_OUT                  ${mol}
ANALYSIS_LMP_TSTEP            ${tStep}
ANALYSIS_MSD_CORLENGTH        0.5
ANALYSIS_MSD_MEMORYMB         2500 
ANALYSIS_MSD_CMCORRECTION     1
ANALYSIS_DIELECTRIC_EWALD     0.001
${grpString}

DATA
