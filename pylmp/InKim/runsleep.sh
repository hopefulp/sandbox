#PBS -l walltime=480:00:00
#PBS -l nodes=1:ppn=12
#PBS -N serial

#!/bin/tcsh

cat $PBS_NODEFILE
set NN = `cat $PBS_NODEFILE | wc -l`

setenv WDIR /qcfs/noische/research/

cd $WDIR
sleep 1000000
#cd /qcfs/noische/research/graphene/solv/kwac/mos2/11p0/NVT
#python3.5 /qcfs/noische/scripts/graphene_mos2/gethbond_parallel.py -b ../*bond.bgf -t *NVT*trj -s "25.21480 <= atom.y <= 87.65070" -o sdmos_gr11p0.nvt
