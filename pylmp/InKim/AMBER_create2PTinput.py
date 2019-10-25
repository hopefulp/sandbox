#!/home/noische/program/python27/bin/python
"""
create2PTinput.py
Original: Dec 28 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import time
import getopt
import re

# Custom Modules
sys.path.append("/home/noische/scripts")
sys.path.append("/home/noische/script")
import nutils as nu

# Globals
version = '120904'

def create2PTinput(crd_file, vel_file, top_file, grps_file, mden_file, ctrl_file, prefix):
	"""
def create2PTinput():
	Write what this function does.

	"""

	### initialize

	### find md_temp and md_avgeng
	l_md_temp = []; l_md_avgeng = [];
	f_mden_file = open(mden_file)
	avg_info = "";
	flag_avg_info = False;
	while 1:
		line = f_mden_file.readline()
		parse = re.split('\s*', line)

		if not line:
			break;

		# average energy
		if 'L0' in parse:
			if "Etot" in parse[3]:
				pass;
			else:
				l_md_avgeng.append(float(parse[3]))
	
		elif 'L1' in line:
			if 'Temp' in parse[1]:
				pass;
			else:
				l_md_temp.append(float(parse[1]))

	# average
	md_temp = reduce(lambda x, y: x + y, l_md_temp) / float(len(l_md_temp))
	md_avgeng = reduce(lambda x, y: x + y, l_md_avgeng) / float(len(l_md_avgeng))


	# old version of finding average temp and Etot
	#avg_info = avg_info.replace("\n", " ")
	#avg_info = avg_info.replace("=", " ")
	#parse = avg_info.split()
	#index_md_temp = parse.index("TEMP(K)") + 1
	#index_md_avgeng = parse.index("Etot") + 1
	#md_temp = float(parse[index_md_temp])
	#md_avgeng = float(parse[index_md_avgeng]) #* 4.1868

	### write ctrl_file
	f_ctrl_file = open(ctrl_file, 'w')
	output = '';

	output += "#structure input options\n"
	output += "IN_AMBERPRMTOP\t\t\t" + str(top_file) + "\n"
	output += "IN_GROUPFILE\t\t\t" + str(grps_file) + "\n"
	output += "\n"

	output += "#2pt options\n"
	output += "ANALYSIS_FRAME_INITIAL\t\t1\n"
	output += "ANALYSIS_FRAME_FINAL\t\t0\n"
	output += "ANALYSIS_FRAME_STEP\t\t1\n"
	output += "ANALYSIS_VAC_CORLENGTH\t\t0.5\n"
	output += "ANALYSIS_VAC_MEMORYMB\t\t10500\n"
	output += "ANALYSIS_VAC_FIXED_DF\t\tg\n"
	output += "ANALYSIS_VAC_LINEAR_MOL\t\t0\n"
	output += "ANALYSIS_VAC_ROTN_SYMMETRY\t2\n"
	output += "ANALYSIS_OUT\t\t\t" + str(prefix) + "\n"
	output += "\n"

	output += "#md options\n"
	output += "MD_AVGTEMPERATURE\t\t" + str(md_temp) + "\n"
	output += "MD_AVGENERGY\t\t\t" + str(md_avgeng) + "\n"
	output += "MD_TSTEP\t\t\t0.002\n"
	output += "\n"

	output += "#trajectory options\n"
	output += "IN_AMBERTRJ\t\t\t" + str(vel_file) + "\n"
	output += "IN_AMBER_COORD_FILE\t\t" + str(crd_file) + "\n"
	output += "TRAJ_DUMPFREQ\t\t\t2\n"

	f_ctrl_file.write(output)
	f_ctrl_file.close()


	### end of create2PTinput


if __name__ == '__main__':

	option = ""; args = ""; crd_file = ""; vel_file = ""; top_file = ""; grps_file = ""; mden_file = ""; ctrl_file = ""; prefix = "";
	number = 0
	usage = """
Usage: create2PTinput.py -c crdfile -v velfile -t topfile -g grpsfile -e mdenfile -i 2pt_ctrlfile -p prefix_4_2pt

	"""

	if len(sys.argv) < 2:
		print(usage); sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hc:v:t:g:e:i:p:', ['help','crd=','vel=','top=','grps=','mden=','ctrl=','prefix='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-c', '--crd'):
			crd_file = value
		elif option in ('-v', '--vel'):
			vel_file = value
		elif option in ('-t', '--top'):
			top_file = value
		elif option in ('-g', '--grps'):
			grps_file = value
		elif option in ('-e', '--mden'):
			mden_file = value
		elif option in ('-i', '--ctrl'):
			ctrl_file = value
		elif option in ('-p', '--prefix'):
			prefix = value
		elif option in (''):
			print(usage); sys.exit(0)

	# default settings

	create2PTinput(crd_file, vel_file, top_file, grps_file, mden_file, ctrl_file, prefix)

