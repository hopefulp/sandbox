#!/home/noische/python

import sys
import string
import copy

import bgf
import bgftools
import nutils as nu

version = 'kdft_20140203'
usage = """
BGF_reorderWaterMolecules.py bgf_name
"""

if len(sys.argv) < 2:
	print(usage)
	sys.exit(0)

myBGF = bgf.BgfFile(sys.argv[1])
newBGF = bgf.BgfFile()

# find atoms for water
aNo_total = set([])
for atom in myBGF.a:
	aNo_total.add(atom.aNo)

aNo_water = set(bgftools.listWaterAtoms(myBGF))

# pick one water molecule
water_O_ffType = ""; water_H_ffType = "";
for ano in aNo_water:
	atom = myBGF.getAtom(ano)
	if "O" in atom.ffType:
		if "W" in atom.ffType or "w" in atom.ffType:
			water_O_ffType = atom.ffType
			atom2 = myBGF.getAtom(atom.CONECT[0])
			if "H" in atom2.ffType:
				if "W" in atom2.ffType or "w" in atom2.ffType:
					water_H_ffType = atom2.ffType

if water_O_ffType != "" and water_H_ffType != "":
	print("Found water FF types: %s %s", (water_O_ffType, water_H_ffType))
else:
	nu.warn("Failed to find water FF type. The script will set FF types as OW and HW for water.")
	water_O_ffType = "OW"
	water_H_ffType = "HW"

# make a new bgf which contains only water
for atom in myBGF.a:
	if atom.aNo in aNo_water:
		atom2 = copy.deepcopy(atom)
		newBGF.addAtom(atom2)
newBGF.renumber()
#newBGF.saveBGF("temp_new.bgf")

# list up indices for water molecules in myBGF
l_delatoms = []
for i in aNo_water:
	l_delatoms.append(myBGF.a2i[i])

# delete water from myBGF
myBGF.delAtoms(l_delatoms)
#myBGF.delAtoms(list(aNo_water))

# renumber myBGF
myBGF.renumber()
#myBGF.saveBGF("temp_my.bgf")

# for atoms in water, write. REMARK: HETATM->ATOM 
for atom in newBGF.a:
	atom.aTag = 0	# change into ATOM
	atom.rName = "WAT"

	if "O" in atom.ffType:
		atom.aName = "O"
		atom.ffType = water_O_ffType

		l_Hatoms = atom.CONECT
		for index, ano2 in enumerate(l_Hatoms):
			atom2 = newBGF.getAtom(ano2)
			atom2.aName = "H" + str(index + 1)
			atom2.ffType = water_H_ffType

# merge myBGF and newBGF
myBGF = myBGF.merge(newBGF, True)
myBGF.renumber()

# save
myBGF = bgftools.renumberMolecules(myBGF, 0)
myBGF.renumber()

try:
	bgf_file = sys.argv[2]
except:
	myBGF.saveBGF(sys.argv[1].split(".bgf")[0] + "_mod.bgf")
else:
	myBGF.saveBGF(bgf_file)

### end
