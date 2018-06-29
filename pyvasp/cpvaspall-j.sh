#!/bin/tcsh
if ($#argv == 0) then
	set nds = 1
	set pwdir = `pwd`
	set title = `basename $pwdir`
else if ($#argv == 1) then
	set nds = $argv[1]
	set pwdir = `pwd`
	set title = `basename $pwdir`
else if ($#argv == 2) then
	set nds = $argv[1]
	set title = $argv[2]
endif

set HOST = `hostname`
switch(${HOST})
        case idft:
                set NPN = 8
                set WTM = 480
                set DQ = dque
        breaksw
        case jdft:
                set NPN = 4
                set WTM = 480
                set DQ = dque
        breaksw
        case kdft:
                set NPN = 12
                set WTM = 480
                set DQ = dque
        breaksw
        case psi:
                set NPN = 16
                set WTM = 120
                set DQ = dque
        breaksw
        case rho:
                set NPN = 16
                set WTM = 120
                set DQ = dque
        breaksw
        case mu:
                set NPN = 32
                set WTM = 120
                set DQ = small
        breaksw
        case qch:
                set NPN = 8
                set WTM = 480
                set DQ = dque
        breaksw
        default:
                set NPN = 8
                set WTM = 480
                set DQ = dque
        breaksw
endsw

if (-e vsp_all.pbs) then
	rm vsp_all.pbs
	cp /qcfs/biduri/scripts/pbs/vsp_all.pbs .
else
	cp /qcfs/biduri/scripts/pbs/vsp_all.pbs .
endif

sed -i -e 's|nodes=N|nodes='${nds}'|' vsp_all.pbs
sed -i -e 's|ppn=NN|ppn='${NPN}'|' vsp_all.pbs
sed -i -e 's|walltime=TTT|walltime='${WTM}'|' vsp_all.pbs
sed -i -e 's|TITLE|'${title}'|' vsp_all.pbs
sed -i -e 's|QQQQ|'${DQ}'|' vsp_all.pbs

sed -i -e 's|TYPE|std|' vsp_all.pbs

if ($HOST == qch) then
	sed -i -e 's|dque|batch|' vsp_all.pbs
endif

if ((`grep 'rlx-xy' INCAR|wc -l`)) then
	sed -i -e 's|vasp-|rlxxy-vasp-|' vsp_all.pbs
endif

if ((`grep 'rlx-orth' INCAR|wc -l`)) then
	sed -i -e 's|vasp-|rlxorth-vasp-|' vsp_all.pbs
endif

if ((`grep 'rlx-z' INCAR|wc -l`)) then
	sed -i -e 's|vasp-|rlxz-vasp-|' vsp_all.pbs
endif

if ((`grep 'k-gamma' INCAR|wc -l`)) then
	sed -i -e 's|-std|-gam|' vsp_all.pbs
endif

if ((`grep 'k-full' INCAR|wc -l`)) then
	sed -i -e 's|-std|-ncl|' vsp_all.pbs
endif
