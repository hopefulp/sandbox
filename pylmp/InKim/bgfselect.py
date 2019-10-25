#!/home/noische/program/python27/bin/python
"""
bgfselect.py
Original: Sep 02 2011 In Kim

Module containing BGF-file selection
"""

# python modules
import sys
import os
import string
import copy
import math

# custom modules
import bgf
import bgftools
import dreiding

version = '110819'

def selectResName(bgf_file, resName, out_file, silent=False):
	"""
selectResName(bgf_file, resName, out_file, silent=False):
	return a BgfFile class which contains the same residue name as designated.
	"""

	# open bgf
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = copy.deepcopy(bgf_file)
	else:
		if not silent: print("Reading " + bgf_file + " ..")
		myBGF = bgf.BgfFile(bgf_file)

	l_notSelected_aNo = []

	# delete unwanted atoms
	for atom in myBGF.a:
		if not resName in atom.rName:
			l_notSelected_aNo.append(atom.aNo)

	l_notSelected_aNo.sort()
	l_notSelected_aNo.reverse()
	myBGF.delAtoms(l_notSelected_aNo)
	if not silent: print(str(len(myBGF.a)) + " atoms are selected.")

	# save
	if isinstance(out_file, str):
		if not silent: print("Saving information to " + out_file + " ..")
		myBGF.saveBGF(out_file)
		return 1;
	else:
		return myBGF;

	### end of selectResName

def selectResNo(bgf_file, resNo, out_file=0, silent=False):
	"""
selectResName(bgf_file, resNo, out_file, silent=False):
	return a BgfFile class which contains the same residue number as designated.
	"""

	resNo = int(resNo)

	# open bgf
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = copy.deepcopy(bgf_file)
	else:
		if not silent: print("Reading " + bgf_file + " ..")
		myBGF = bgf.BgfFile(bgf_file)

	l_notSelected_aNo = []

	if not silent: print("entering selection..")

	# delete unwanted atoms
	for atom in myBGF.a:
		if not resNo == atom.rNo:
			l_notSelected_aNo.append(myBGF.a2i[atom.aNo])

	#l_notSelected_aNo.sort()
	#l_notSelected_aNo.reverse()
	myBGF.delAtoms(l_notSelected_aNo)
	if not silent: print(str(len(myBGF.a)) + " atoms are selected.")

	# save
	if isinstance(out_file, str):
		if not silent: print("Saving information to " + out_file + " ..")
		myBGF.saveBGF(out_file)
		return 1;
	else:
		return myBGF;

	### end of selectResName


#-------------------------------------
#
# - Print information
#
if __name__ == '__main__':

	# get directory:
	directory = dir()

	# set imported stuff we don't want to see:
	imported = ['sys', 'os', 'string', 'copy', 'math', 'types', 'bgf', 'bgftools', 'dreiding']

	# print __doc__ for the module:
	print("\n")
	print("-"*60)
	if 'version' not in directory:  version = '??????'
	print("%-45s%15s" % (os.path.basename(sys.argv[0]), 'ver: '+version))

	print("-"*60)
	print(__doc__)


	# import types:
	import types

	# create hash-table:
	hashtable = {}
	for item in directory:
		actual_item = eval(item)
		if item in imported:
			# don't show imported stuff:
			pass
		elif type(actual_item) is types.ModuleType:
			# don't discuss other modules:
			pass
		elif type(actual_item) is types.FunctionType:
			# show __doc__s for functions:
			hashtable[item] = actual_item.__doc__
		elif type(actual_item) is types.ClassType:
			# show __doc__s for classes:
			title = item+' class: '
			hashtable[item] = title +  ( '-' * (60-len(title)) )
			hashtable[item] += actual_item.__doc__

			# show __doc__s for class elements:
			for classItem in dir(actual_item):
				actual_class_item = eval(item+'.'+classItem)
				if type(actual_class_item) is types.ModuleType:
					# don't discuss other modules:
					pass
				elif type(actual_class_item) is types.UnboundMethodType \
						or type(actual_class_item) is types.MethodType:
					# show __doc__s for functions:
					hashtable[item] += actual_class_item.__doc__ + '\n'

				elif classItem in ['__doc__','__module__']:
					pass
				else:
					# for other stuff show the value:
					hashtable[item] += '\n'+classItem+' = '+str(actual_class_item)+'\n'

			hashtable[item] +=  ( '-'*60 )+'\n\n'
			
		elif item[0] != '_':
			# for other stuff show the value:
			hashtable[item] = '\n'+item+' = '+str(actual_item)+'\n'

	# print info out
	keys = hashtable.keys()
	keys.sort()
	
	print("Contents:")
	print("-"*60 )
	for item in keys:
		print(hashtable[item])

	print("\n")
	print("-"*60 )
	print("contact: noische@kaist.ac.kr")
	print("-"*60 )
	print("\n")
	
	# done!
