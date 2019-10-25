#!/home/noische/python

import sys
import os
import os.path
import bgf
import bgftools
from PEI_isomorphism_serial import *

def generatePolymer(n_generation):
	"""
	"""

	curr_dir = os.path.abspath(".") + "/"
	products_dir = curr_dir + "generated/"
	if not os.path.exists(products_dir):
		os.makedirs(products_dir)

	#generate_command = '/home/noische/scripts/buildPolymer.pl -m "primary.bgf secondary.bgf tertiary.bgf" -f /home/noische/ff/dreiding-den.par -w 600 -r "35 35 30" -s '
	generate_command = '/home/noische/scripts/buildPolymer.pl -m "/qcfs/noische/research/PEI/PEI-600/b3lyp/monomer_b3lyp_20120615/primary.b3lyp.bgf /qcfs/noische/research/PEI/PEI-600/b3lyp/monomer_b3lyp_20120615/secondary.b3lyp.final.bgf /qcfs/noische/research/PEI/PEI-600/b3lyp/monomer_b3lyp_20120615/tertiary.b3lyp.final.bgf" -f /home/noische/ff/dreiding-den.par -w 600 -r "35 35 30" -s '

	n_count = 0;	# increase when success
	n_try = 1;	# number of trials

	while n_count <= n_generation:

		bgffile = "600_" + "{:0=4}".format(n_count) + ".bgf"
		if os.path.isfile(products_dir + bgffile):
			n_count += 1;
			continue;

		run_command = generate_command + bgffile
		#run_command += "> /dev/null"	# postfix

		print("\n\n***** Building a polymer.. (success " + str(n_count) + " / try " + str(n_try) + ")")
		print("")
		os.system(run_command)

		# check amine ratio
		aminegroup = bgftools.getAmineGroupInfo(bgffile)
		print(" Number of amine groups in the generated polymer: " + str(aminegroup))
		if aminegroup == [5, 5, 4]:
			print("\t** Passed amine group ratio test. Moving to isomorphism test.")
		else:
			print("\t** Build Failed. (" + str(aminegroup) + "). File will be deleted and regenerated.")
			os.system("rm " + bgffile)
			os.system("rm -rf files")
			n_try += 1;
			continue;

		# check isomorphism
		iso = check_isomorphism(products_dir, bgffile, True)	# check isomorphism of all structures
		if not iso == True:
			print("\t** Build Failed. (the generated structure is identical to " + iso + " ) File will be deleted and regenerated.")
			os.system("rm " + bgffile)
			os.system("rm -rf files")
			n_try += 1;
			continue;

		# move the generated file to ./generated
		src = curr_dir + bgffile
		dst = products_dir + bgffile
		os.rename(src, dst)
		print("***** Build successful.")
		print("")
		n_count += 1;

		os.system("rm -rf files")
		n_try += 1;

		# end of while

if __name__ == "__main__":
	n = int(sys.argv[1])
	generatePolymer(n)
