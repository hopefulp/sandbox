#!/bin/bash
## usage $0 job dir version

qcvi=/qcfs/joonho/VaspINI
qcvin=/qcfs/joonho/binvasp

mag=$1
metal=$2
cont=$3
prec=$4
gga=$5
opt=$6
log=$7

if [ $# -ne 7 ]; then
    echo "Usage:: $0 mag metal cont prec gga opt log"
    echo "  mag:: AFM FM NM"
    echo " cont:: ini wav chg"
    echo " prec:: prec"
    echo "  gga:: rpd rpdu"
    echo "  opt:: sp opt crelax"
    echo "  log:: log for default, logall for all electron writing"
    exit
fi

$qcvin/vrun_incar.pl $qcvi/INC/INCAR.1.$mag $metal

cat  INCAR.mag > INCAR
cat  $qcvi/INC/INCAR.2.$cont >> INCAR
cat  $qcvi/INC/INCAR.3.$prec >> INCAR
cat  $qcvi/INC/INCAR.4.$gga  >> INCAR
cat  $qcvi/INC/INCAR.5.$opt >> INCAR
cat  $qcvi/INC/INCAR.6.$log >> INCAR


