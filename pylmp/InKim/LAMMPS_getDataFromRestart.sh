#!/bin/bash

# this script reads LAMMPS restart file and writes coordinates to data file
# usage: LAMMPS_getDataFromRestart.sh in_singlepoint_file restart_file output_data_file working_directory
# $1 in_singlepoint_file
# $2 restart_file
# $3 output_data_file
# $4 working directory

usage="usage: LAMMPS_getDataFromRestart.sh in_singlepoint_file restart_file output_data_file working_directory"
if [ $# -lt 3 ]
then
	echo $usage
    exit
fi

recipe="/qcfs/noische/scripts/dat/LAMMPS/in.lammps.write_data_from_restart"
workdir=$4
workfile=$workdir/_in.getData

echo "LAMMPS_getDataFromRestart.sh: cd $4"
cd $4
cat $1 $recipe > $workfile
sed -i "s/run 0//" $workfile
sed -i -r "s/read_data.*/read_restart $2/" ${workfile}
sed -i -r "s/write_data.*/write_data $3/" ${workfile}
echo "Running LAMMPS to fetch coordinates from the restart file: $EXEC"
$EXEC -in $workfile -screen none
#rm _in.getData # log.cite log.lammps
echo "Done. Check $3"
