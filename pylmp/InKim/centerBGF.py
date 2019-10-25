#!/home/noische/program/python27/bin/python

import sys
import re
import string
import getopt
import optparse
from os import popen

# BGF modules
sys.path.append("/home/noische/script")
import bgf
import bgftools
import nutils as nu

#-----------------
# update the coordinate in the original BGF file from LAMMPS trajectory file
#_________________
def centerbgf(bgf_file, out_file, ff_file, method = "com_origin", silent=False):
	if silent == 1:
		nu.shutup()

	boxcenter = [];

	# open bgf
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print("Reading " + bgf_file + " ..")
		myBGF = bgf.BgfFile(bgf_file)
		if not silent: print("Moving the origin of the coordinate in " + bgf_file + " and saving to " + out_file + ".")

	if method != "box_origin":
		new_x, new_y, new_z = bgftools.getCom(myBGF, ff_file)

	if len(myBGF.CRYSTX) > 2:
		boxcenter = [ (i/2.0) for i in myBGF.CRYSTX[0:3] ]
	else:
		boxsize = bgf.getBGFSize(myBGF, 0)		# [xlo, xhi, ylo, yhi, zlo, zhi]
		boxcenter = [ (boxsize[1] - boxsize[0])/2, (boxsize[3] - boxsize[2])/2, (boxsize[5] - boxsize[4])/2 ]

	if method == "com_origin":
		bgf.moveBGF(myBGF, (-new_x), (-new_y), (-new_z))
	elif method == "com_center":
		dx = boxcenter[0] - new_x
		dy = boxcenter[1] - new_y
		dz = boxcenter[2] - new_z
		bgf.moveBGF(myBGF, dx, dy, dz)
	elif method == "box_origin":
		dx = boxcenter[0]
		dy = boxcenter[1]
		dz = boxcenter[2]
		bgf.moveBGF(myBGF, dx, dy, dz)
	else:
		bgf.moveBGF(myBGF, (-new_x), (-new_y), (-new_z))

	# save
	if isinstance(out_file, str):
		if not silent: print("Saving information to " + out_file + " ..")
		myBGF.saveBGF(out_file)
		return 1;
	else:
		return myBGF;


def centerbgffileclass(myBGF, method = "com_origin"):

	new_x, new_y, new_z = bgftools.getCom(myBGF, "")

	if len(myBGF.CRYSTX) > 2:
		boxcenter = [ (i/2.0) for i in myBGF.CRYSTX[0:3] ]
	else:
		boxsize = bgf.getBGFSize(myBGF, 0)
		boxcenter = [ (boxsize[1] - boxsize[0])/2, (boxsize[3] - boxsize[2])/2, (boxsize[5] - boxsize[4])/2 ]

	if method == "com_origin":
		bgf.moveBGF(myBGF, (-new_x), (-new_y), (-new_z))
	elif method == "com_center":
		dx = boxcenter[0] - new_x
		dy = boxcenter[1] - new_y
		dz = boxcenter[2] - new_z
		bgf.moveBGF(myBGF, dx, dy, dz)
	elif method == "box_origin":
		dx = boxcenter[0]
		dy = boxcenter[1]
		dz = boxcenter[2]
		bgf.moveBGF(myBGF, dx, dy, dz)
	else:
		bgf.moveBGF(myBGF, (-new_x), (-new_y), (-new_z))

	return myBGF;


if __name__ == "__main__":

	option = ""; args = ""; bgf_file = ""; ff_file = ""; out_file = ""; method = "";
	usage = """
centerBGF: read coordinates data from the BGF file
           move the coordinates as desired method
           write the data to the targeted BGF file
Usage: centerBGF.py -b bgf_file -t trj_file -o out_file
	"""

	options, args = getopt.getopt(sys.argv[1:], 'hb:f:o:m:', ['help','bgf=','ff=','out=','method='])
	print("Requested options: " + str(options))
	for option, value in options:
		if option in ('-h', '--help'):
			print usage
			sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-f', '--ff'):
			ff_file = value
		elif option in ('-o', '--out'):
			out_file = value
		elif option in ('-m', '--method'):
			method = value
		elif option == NULL:
			print usage
			sys.exit(0)

	# main call
	centerbgf(bgf_file, ff_file, out_file, method)

