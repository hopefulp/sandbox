#!/home/noische/program/python27/bin/python

# System Modules
import sys
import re
import string
import getopt
import optparse
import math
import time
from os import popen

# BGF modules
sys.path.append("/home/noische/script")
import bgf
import nutils as nu
import bgftools

# Pizza.py modules
sys.path.append("/home/noische/program/pizza-1Oct10")
import dump

def sortkey(list):
	return list[0];

#-----------------
# update the coordinate in the original BGF file from LAMMPS trajectory file
#_________________
def updatebgf(bgf_file, trj_file, out_file, periodic, silent):
	if not silent: print("Updating " + bgf_file + " with the coordinate information in " + trj_file + " to " + out_file + ".")

	# load LAMMPS trajectory
	myDUMP = dump.dump(trj_file)
	#myDUMP.unwrap()

	l_timestep = [];
	# for bypassing irrelevant trajectories
	#	while 1:
	#		time_myDUMP = myDUMP.next()
	#		if not silent: sys.stdout.write("\rBypassing the trajectory " + str(time_myDUMP))
	#		if time_myDUMP == -1:
	#			break;
	#		l_timestep = myDUMP.time()
	#
	if not silent: sys.stdout.write("\n")

	l_timestep = myDUMP.time()
	if len(l_timestep) == 0:
		nu.die("No trajectory found on " + trj_file + " !!")
	elif len(l_timestep) >= 2:
		myDUMP.tselect.one(l_timestep[-1])
	else:
		myDUMP.tselect.one(l_timestep[0])

	myDUMP.sort()
	time, box, atoms, bonds, tris, lines = myDUMP.viz(0)
	atoms.sort(key=sortkey)
	atomcoord = dict();

	for item in atoms:
		atomcoord[str(int(item[0]))] = [item[2], item[3], item[4]]

	# load BGF
	myBGF = bgf.BgfFile(bgf_file)

	# Update coordinate information in BGF File
	if len(myBGF.a) != len(atoms):
		nu.die("The number of atoms in " + bgf_file + " and " + trj_file + " is different.")
	else:
		# update the coordinate information
		for myatom in myBGF.a:	# using atoms.sort(key=sortkey)
			myatom.x = atomcoord[str(myatom.aNo)][0]
			myatom.y = atomcoord[str(myatom.aNo)][1]
			myatom.z = atomcoord[str(myatom.aNo)][2]

	# Update periodic information in BGF File if exists
	if myBGF.PERIOD != "" and myBGF.CRYSTX != []:
		if not silent: print("Writing periodic box information..")
		#print(str(box))
		myBGF.CRYSTX[0] = float(box[3]) - float(box[0])
		myBGF.CRYSTX[1] = float(box[4]) - float(box[1])
		myBGF.CRYSTX[2] = float(box[5]) - float(box[2])

	# apply periodic condition if periodic == True
	#if periodic:
	#	myBGF = bgftools.periodicMoleculeSort(myBGF, 0)

	myBGF.saveBGF(out_file)

	if not silent: print("The BGF file " + bgf_file + " is updated with the coordinates in " + trj_file + ".")

	return 1;	# success

	### end of updatebgf


if __name__ == "__main__":

	option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; step = 0;
	usage = """
	updateBGF: read coordinate data from the LAMMPS trajectory file
	           write the data to the original BGF file
	Usage: updateBGF.py -b bgf_file -t trj_file -o out_file
	"""
	
	periodic = False;

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(1)

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:o:p', ['help','bgf=','trj=','out='])
	print("Requested options: " + str(options))
	for option, value in options:
		if option in ('-h', '--help'):
			print(usage)
			sys.exit(0);
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-t', '--trj'):
			trj_file = value
		elif option in ('-p', '--periodic'):
			periodic = True
		elif option in ('-o', '--out'):
			out_file = value
		elif option == NULL:
			print(usage)
			sys.exit(0)

	# main call
	updatebgf(bgf_file, trj_file, out_file, periodic, silent=False)

