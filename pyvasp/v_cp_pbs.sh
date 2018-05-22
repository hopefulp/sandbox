#!/bin/bash
### modified from cpvaspall.sh 
if [ $# -eq 0 ]; then
    echo  Input job name for directory. So, exit.
    exit
elif [ $# -eq 1 ]; then
	title=$1
    cvasp="std"
	nds=1
elif [ $# -eq 2 ]; then
	title=$1
    cvasp=$2
	nds=1
fi
echo jobname \"$title\" :  vasp-type \"$cvasp\"
HOST=`hostname`
case ${HOST} in
    "idft")
        NPN=8
        WTM=480
        DQ=dque
        ;;
    "jdft")
        NPN=4
        WTM=480
        DQ=dque
        ;;
    "kdft")
        NPN=12
        WTM=480
        DQ=dque
        ;;
    "psi")
        NPN=16
        WTM=120
        DQ=dque
        ;;
    "rho")
        NPN=16
        WTM=120
        DQ=dque
        ;;
    "mu")
        NPN=32
        WTM=120
        DQ=small
        ;;
    "qch")
        NPN=8
        WTM=480
        DQ=dque
        ;;
    *)
        NPN=8
        WTM=480
        DQ=dque
        ;;
esac

echo "cp from /qcfs/joonho/bin/template/vsp_all.pbs"
cp /qcfs/joonho/bin/template/vsp_all.pbs .

sed -i -e 's|nodes=N|nodes='${nds}'|' vsp_all.pbs
sed -i -e 's|ppn=NN|ppn='${NPN}'|' vsp_all.pbs
sed -i -e 's|walltime=TTT|walltime='${WTM}'|' vsp_all.pbs
sed -i -e 's|TITLE|'${title}'|' vsp_all.pbs
sed -i -e 's|QQQQ|'${DQ}'|' vsp_all.pbs

sed -i -e 's|TYPE|'${cvasp}'|' vsp_all.pbs

