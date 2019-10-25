#!/home/noische/python
"""
clusterSetting.py
Original: Aug 10 2011 In Kim
"""

# Python Modules
import sys
import os

# Globals
version = "120625"
cluster = "kdft"
usage = """
clusterSetting.py: Do not access this file from the command line.
"""

scratch_dir = "/scratch/noische"
hostname = os.environ['HOSTNAME']
n_proc = 0;

try:
	hostnum = int(hostname.split('kdft')[1])
except:
	n_proc = 1;	# headnode
else:
	if hostnum < 25:
		n_proc = 12;
	else:
		n_proc = 16;

#mpi_command = "/opt/mpi/intel/mpich2-1.4rc2/bin/mpirun -f hostfile -np " + str(n_proc)
#lammps_command = "/opt/applic/lammps/bin/lmp_kdft "
mpi_command = "opt/intel/composer/composer_xe_2013.1.117/mpi/openmpi-1.6.3/bin/mpirun -hostfile $PBS_NODEFILE -n $NPROC " + str(n_proc)
lammps_command = "/qcfs/noische/program/kdft/bin/lmp_kdft "

#file = open("hostfile", 'w')
#for i in range(n_proc):
#	file.write(os.environ['HOSTNAME'] + "\n")
#file.close()

### end of loadClusterSetting

if __name__ == "__main__":

	print(usage)
	sys.exit(0)
