#!/bin/bash

## Usage:: V_run_ini.sh   ac[ea|ey|pa|py|zz] Sc AFM[FM|nonmag] KPOINT  af[fm|nm]

if [ $# -ne 5 ]; then
    echo "Usage:: $0  POS metal INC KPOINT dir-suffix"
    exit
fi

if [ $1 == "usage" -o $1 == "Usage"]; then
    echo "     :: POS from POSCAR.Fe.pos where pos is \"zz, ac, py, pa etc\" "
    echo "     :: INC from incar.AFM or FM or nonmag"
    echo "     :: KPOINTS is the number for k-points"
    echo "     :: dir-suffix is added as magnetism of af, fm, nonmag"
    exit
fi

qcvin=/qcfs/joonho/binvasp
qcvi=/qcfs/joonho/VaspINI

pos=$1			# position from zz ac 
metal=$2		# 
mag=$3			# change magnetism, FM, AFM, AllFM
kpoints=$4		# gamma point
dirsuff=$5

dir=$metal-$pos-$dirsuff

if [ $metal == "Sc" -o $metal == "Ca" ]; then
    pot_id=${metal}sv
else
    pot_id=${metal}
fi


mkdir $dir
### POSCAR
$qcvin/vrun_poscar.pl $qcvi/POSCAR.Fe.$pos $metal
mv POSCAR $dir
### POTCAR
cp $qcvi/Pot/pawpbe_${pot_id}OCH.pot $dir/POTCAR
### INCAR			$metal for magnitude of magnetization
$qcvin/vrun_incar.sh $mag $metal ini prec rpd crelax log
cp INCAR $dir

### KPOINTS
cp $qcvi/kp1.gamma	$dir/KPOINTS

$qcvin/changeline_pbs_$HOST.pl pbs-$HOST.csh $dir



