#!/home/noische/program/python27/bin/python
"""
example.py
Contains some examples dealing with bgffile and bgfatom classes.
"""

import sys
if '/home/noische/script' not in sys.path:
	sys.path.append('/home/noische/script')
from bgf import *

myBGF = BgfFile(sys.argv[1])
#print("Opened BGF file:")
#print(myBGF)
#print("\n")

# declare a new atom, 'newatom'
newatom = BgfAtom()

# change properties of 'newatom'
newatom.aTag = 1
newatom.aNo = 11
newatom.aName = "C11  "
newatom.rName = "PRI"
newatom.chain = "X"
newatom.rNo = 0
newatom.x = 1.0
newatom.y = 2.0
newatom.z = 3.0
newatom.ffType = "C_3  "
newatom.bonds = 4
newatom.lpair = 0
newatom.charge = 0.00000

# you can print newatom by print() in a bgf format.
#print("Newly added atom: newatom")
#print(newatom)
#print("\n")

# add the new atom to first
myBGF.addAtom(newatom)

# internal index can be accessed using a2i dictionary.
# if the function requires an index, you have to change the atom number by calling a2i[].
#print("Interal index of newatom:")
#print(myBGF.a2i[newatom.aNo])

# check if two atoms are bonded?
# myBGF.a is a list that contains atoms
#print("Are they bonded? " + str(is_bonded(newatom, myBGF.a[0])))

# connect two atoms by throwing two atom Numbers to connectAtoms()
# the following function is same as 
# 1) myBGF.connect(myBGF.a2i[1], myBGF.a2i[11])
# 2) myBGF.connect(myBGF.getAtomIndex[1], myBGF.getAtomIndex[11])
myBGF.connectAtoms(1, 11)
print(myBGF.getAtom(15))

# remove bonds between two atoms by disconnect()
myBGF.disconnect(myBGF.a2i[1], myBGF.a2i[11])

# changed BGF file contains this information
#print("*"*10 + " Changed BGF file: " + "*"*10)
#print(myBGF)
#print("\n")

# save the changed myBGF structure
#myBGF.saveBGF(sys.argv[1] + "_new")
exit
