#!/home/noische/python

import sys
import os
import string
import copy
import time

import bgf
import bgftools
import nutils as nu

usage = """
PEGDE_generator.py head head-monomer monomer tail-monomer tail target_mw ff_file

in head bgf, the following should be done:
	- H_EG is the hydrogen atom that H_H in monomer atom will be located.
	- The chain of head C which briges the monomer and head should be "H".
in monomer bgf, the following should be done:
	- H_H is the hydrogen attached at the head Carbon of chain H 
	- H_T is the hydrogen attached at the tail Oxygen of chain T
in tail bgf, the following should be done:
	- H_EG is the hydrogen atom that H_T in monomer atom will be located.
	- The chain of head C should be "T".
in ff_file, you can use "" to load default ff, DREIDING 2.

"""


def main():

	if len(sys.argv) < 2:
		print usage;
		sys.exit(0)
	
	head = bgf.BgfFile(sys.argv[1]);
	head_mon = bgf.BgfFile(sys.argv[2]);
	monomer = bgf.BgfFile(sys.argv[3]);
	tail_mon = bgf.BgfFile(sys.argv[4]);
	tail = bgf.BgfFile(sys.argv[5]);
	target_mw = float(sys.argv[6]);
	ff_file = sys.argv[7];

	# number of required monomers to make a PEGDE chain
	mw_head = bgftools.getMoleculeMass(head, ff_file)
	print("** Mw of Molecule head: " + "{0:5.3f}".format(mw_head))

	mw_tail = bgftools.getMoleculeMass(tail, ff_file)
	print("** Mw of Molecule tail: " + "{0:5.3f}".format(mw_tail))

	mw_monomer = bgftools.getMoleculeMass(monomer, ff_file)
	print("** Mw of Molecule monomer: " + "{0:5.3f}".format(mw_monomer))

	print("** Requested Mw: " + str(target_mw))

	# number of required monomers
	n_monomers = round((target_mw - mw_head - mw_tail) / mw_monomer)
	print("** Required number of monomers: " + str(int(n_monomers)))

	actual_mw = mw_head + mw_tail + mw_monomer * n_monomers
	print("** Actual Mw: " + "{0:5.3f}".format(actual_mw))


	### overall structure
	### head -- head_mon -- monomer -- tail_mon -- tail


	### connect blocks
	m1 = copy.deepcopy(head_mon)	# start adding monomers from the head_near monomer
	for atom in m1.a:
		atom.rNo = 1;
		atom.rName = 'MON';

	print("* attaching head_monomer")
	n_monomers = int(n_monomers)
	
	### monomer - monomer: 1H_T -- 2H_H
	for i in range(n_monomers - 1):

		# if the last monomer, attach tail_monomer
		# else prepare monomer 2
		if i == n_monomers-2:
			print("* attaching tail_monomer")
			m2 = copy.deepcopy(tail_mon)
		else:
			print("* attaching monomer " + str(i))
			m2 = copy.deepcopy(monomer)
		
		# renumber residue numbers
		for atom in m2.a:
			atom.rNo = i+2
			atom.rName = 'MON';

		# move H_H of m2 to H_T of m1 position
		h_t = 0; h_h = 0;
		for atom in m1.a:
			if atom.ffType == "H_T" and atom.rNo == i+1:
				h_t = atom
				
		for atom in m2.a:
			if atom.ffType == "H_H" and atom.rNo == i+2:
				h_h = atom

		for atom in m2.a:
			atom.x += h_t.x
			atom.y += h_t.y
			atom.z += h_t.z

		# delete H_T in monomer
		m1.delAtom(m1.a2i[h_t.aNo])
		m1.renumber()

		# delete H_H in monomer2
		m2.delAtom(m2.a2i[h_h.aNo])
		m2.renumber()

		# merge
		m1 = m1.merge(m2, True)
		m1.renumber()

		# find tag in the merged molecule
		atomH = bgf.BgfAtom();
		atomT = bgf.BgfAtom();
		for atom in m1.a:
			if 'T' in atom.chain and atom.rNo == i+1:
				atomT = atom
			if 'H' in atom.chain and atom.rNo == i+2:
				atomH = atom

		# update chain
		atomT.chain = "A"; atomH.chain = "A"

		# connect
		m1.connect(m1.a2i[atomH.aNo], m1.a2i[atomT.aNo])


	#m1.saveBGF('temp_monomer.bgf')


	### head - monomer: H_EG -- H_H
	# head
	h1 = copy.deepcopy(head)
	for atom in h1.a:
		atom.rNo = 100;
		atom.rName = 'HED'

	h_eg = 0; h_h = 0;
	for atom in h1.a:
		if atom.ffType == "H_EG":
			h_eg = atom
	
	for atom in m1.a:
		if atom.ffType == "H_H":
			h_h = atom

	for atom in m1.a:
		atom.x += h_eg.x
		atom.y += h_eg.y
		atom.z += h_eg.z

	h1.delAtom(h1.a2i[h_eg.aNo])
	h1.renumber()

	m1.delAtom(m1.a2i[h_h.aNo])
	m1.renumber()

	h1 = h1.merge(m1, True)
	h1.renumber()

	atomH = bgf.BgfAtom(); atomT = bgf.BgfAtom()
	for atom in h1.a:
		if 'T' in atom.chain and 'HED' in atom.rName:
			atomT = atom
		if 'H' in atom.chain and 'MON' in atom.rName:
			atomH = atom

	h1.connect(h1.a2i[atomT.aNo], h1.a2i[atomH.aNo])

	#h1.saveBGF('temp_h-m.bgf')


	# monomer - tail: H_T -- H_EG
	t1 = copy.deepcopy(tail)
	for atom in t1.a:
		atom.rNo = 200;
		atom.rName = 'TAL'

	h_t = 0; h_eg = 0;
	for atom in h1.a:
		if atom.ffType == "H_T":
			h_t = atom

	for atom in t1.a:
		if atom.ffType == "H_EG":
			h_eg = atom

	for atom in t1.a:
		atom.x += h_t.x
		atom.y += h_t.y
		atom.z += h_t.z

	h1.delAtom(h1.a2i[h_t.aNo])
	h1.renumber()

	t1.delAtom(t1.a2i[h_eg.aNo])
	t1.renumber()

	h1 = h1.merge(t1, True)
	h1.renumber()

	atomH = bgf.BgfAtom(); atomT = bgf.BgfAtom()
	for atom in h1.a:
		if 'T' in atom.chain and 'MON' in atom.rName:
			atomT = atom
		if 'H' in atom.chain and 'TAL' in atom.rName:
			atomH = atom

	h1.connect(h1.a2i[atomT.aNo], h1.a2i[atomH.aNo])

	# adding remarks
	h1.REMARK.append("Polymer Generated at " + str(time.asctime(time.gmtime())))
	h1.REMARK.append("head: " + sys.argv[1] + ", head_mon: " + sys.argv[2] + ", monomer: " + sys.argv[3] + ", tail_mon: " + sys.argv[4] + ", tail: " + sys.argv[5])
	h1.REMARK.append("Target Mw: " + str(sys.argv[6]) + ", Actual Mw: " + "{0:5.3f}".format(actual_mw))
	h1.REMARK.append("Number of monomers: " + str(n_monomers))
	h1.REMARK.append("ff: " + sys.argv[7])
	h1.saveBGF('generated_charge.bgf')

	
main()

