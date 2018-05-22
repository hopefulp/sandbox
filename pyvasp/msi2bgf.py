#!/home/jackjack5/epd/bin/python

# 26Apr2017 developed by joonho
from __future__ import print_function
import argparse
import sys
import os
import re
import math
import Bgf_Atom

parser = argparse.ArgumentParser()
explain=" input .msi file -  msi file reads atom and bond with molecular unit "
parser.add_argument("msi_file", help=explain)
args = parser.parse_args()

if re.search("msi$", args.msi_file):
    prefix=args.msi_file[:-4]
    msifile = prefix+'.msi'
    bgffile = prefix+'.bgf'
    print ("convert ", msifile, " to ", bgffile)
else:
    print ("input .msi file")
    exit(0)

inf  = open(msifile)
outf = open(bgffile, 'w')

Remark = 'BGF file created by Joonho'
FF     = 'DREIDING'
Period = '111'
#calculate a,b,c,alpha,beta,gamma from lattice vector
pi = 3.14159265
def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))
def length(v):
	return math.sqrt(dotproduct(v, v))
def angle(v1, v2):
	return (math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))) / pi * 180

while 1:
    # try to find FF, periodicity
    line = inf.readline()
    if not line: break
    if 'A3' in line:
        parse = line.split()
        ax = float(parse[3][1:])
        ay = float(parse[4])
        az = float(parse[5][:-2])
        veca = [ax,ay,az]
    elif 'B3' in line:
        parse = line.split()
        bx = float(parse[3][1:])
        by = float(parse[4])
        bz = float(parse[5][:-2])
        vecb = [bx,by,bz]
    elif 'C3' in line:
        parse = line.split()
        cx = float(parse[3][1:])
        cy = float(parse[4])
        cz = float(parse[5][:-2])
        vecc = [cx,cy,cz]
    elif 'Atom' in line:
        break

la = length(veca)
lb = length(vecb)
lc = length(vecc)
alpha = angle(veca,vecb)
beta  = angle(vecb,vecc)
gamma = angle(vecc,veca)

inf.close()
inf = open(msifile)

str_fprint1 = '%6s %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n'
str_fprint2 = '%-6s %5s %-5.5s %3s %1s %5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f\n'
#### the new format cannot be applied to file.write
#str_print1  = '{:6} {:11.5f}{:11.5f}{:11.5f}{:11.5f}{:11.5f}{:11.5f}'
#str_print2  = '{:6s} {:>5} {:5.5} {:3} {:1} {:^5}{:10.5f}{:10.5f}{:10.5f} {:<5}{:>3d}{:>2d} {:8.5f}'

outf.write('XTLGRF 200 \n')                             # (A6,I5) I5 for version of biogrf
outf.write('DESCRP ' + prefix + '\n')                   # ('DESCRP',1X,A8) short descrp for file
outf.write('REMARK ' + Remark + '\n')                   # ('REMARK',1X,A) File Description - remark,
outf.write('FORCEFIELD ' + FF + '\n')                   # ('FORCEFIELD',1X,A8)
outf.write('PERIOD ' + Period + '\n')                   # ('PERIOD',1X,3I1) 
outf.write('AXES   ZYX\n')                              # ('AXES',3X,A)
outf.write('SGNAME P 1                  1    1\n')      # ('SGNAME',1X,A8,1X,A8,2I5)
#print str_print1.format('CRTSTX',la,lb,lc,alpha,gamma,beta)   #('CRYSTX',1X,6F11.5)
outf.write(str_fprint1 % ('CRTSTX',la,lb,lc,alpha,gamma,beta)) #('CRYSTX',1X,6F11.5)
outf.write('CELLS    -1    1   -1    1   -1    1\n')            # ('CELLS',1X,6I5) ignores in reading
outf.write('FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5, i3,i2,1x,f8.5)\n')

Key_atom='Atom'
Key_bond='Bond'
tag_atom="NO"
tag_bond="NO"
i=0
# to save atoms
hetatm=[]
# for making connection
msi2id={}           # key:value = msi_id:atom_id
bonds=[]
natoms=0
nbonds=0
while 1:
    i+=1
    line = inf.readline()            # read ACL
    if not line: break
    #print line
    if tag_atom == "NO" and tag_bond == 'NO':
        if Key_atom in line:
            #### make dictionary of msi_id:atom_id
            col=re.split('\W+', line)
            msi_id   = col[1]
            natoms  += 1
            tag_atom = 'YES'
            
        elif Key_bond in line:
            nbonds  += 1
            tag_bond = 'YES'
        continue
    elif tag_atom == 'YES':
        if not re.search("\w", line):   # if not word at the end of block : ")"
            tag_atom = 'NO'
            #print str_print2.format('HETATM', atom_index, atom_label, 'RES', 'A','444',x, y, z, atom_fftype, max_cb, num_lp, atom_chg)
            # do not write but save for mixed order in msi file
            #outf.write(str_fprint2 % ('HETATM', atom_index, atom_label, 'RES', 'A',' 444 ',x, y, z, atom_fftype, max_cb, num_lp, atom_chg))
            hetatm_line=[atom_index, atom_label, x, y, z, atom_fftype, max_cb, num_lp, atom_chg]
            hetatm.append(hetatm_line)
            continue
        if 'ACL' in line:
            col = re.split("\W+",line)
            atom_num = col[4]
            atom_sym = col[5]
        elif 'XYZ' in line:
            col = re.split("\s+",line)
            x   = float(col[4].replace('(',""))
            y   = float(col[5])
            z   = float(col[6].replace(')',""))
        elif 'Id' in line:
            col = re.split("\W+",line)
            atom_index = col[4]
            atom_label = atom_sym + atom_index
            msi2id[msi_id]=atom_index
        elif 'Charge' in line:
            col = re.split("\s+",line)
            atom_chg = float(col[4].replace(')',""))
        elif 'FFType' in line:
            col = re.split("\W+",line)
            atom_fftype = col[4]
            max_cb=Bgf_Atom.get_max_CB(atom_fftype)
            num_lp=Bgf_Atom.get_num_lp(atom_fftype)
        continue
    elif tag_bond == 'YES':
        if not re.search("\w", line):
            tag_bond = 'NO'
            continue
        if 'Atom1' in line:
            col = re.split("\W+",line)
            id1=msi2id[col[4]]
        elif 'Atom2' in line:
            col = re.split("\W+",line)
            bonds.append([id1, msi2id[col[4]]])
        continue
#### End of file

#### sort hetatm[] and write
hetatm.sort(key=lambda x:int(x[0]))
#print (*hetatm, sep='\n')
#str_fprint2 = '%-6s %5s %-5.5s %3s %1s %5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f\n'
for fatom in hetatm:
    outf.write(str_fprint2 % ('HETATM',fatom[0],fatom[1],'RES','A',' 444 ',fatom[2],fatom[3],fatom[4],fatom[5],fatom[6],fatom[7],fatom[8]))

#### msi file keeps Atom and Bond in molecular order
outf.write("FORMAT CONECT (a6,12i6)\n")
ai_bonds=[]                 # atom_index is list_index + 1, (id:1) == ai_bonds[0]
for i in range(natoms):     # prepair empty array in advance
    ai_bonds.append([])
#### makes connection w.r.t. atom_index
while bonds:
    [b1s, b2s]=bonds.pop(0)
    b1=int(b1s)
    b2=int(b2s)

    ai_bonds[b1-1].append(b2)   # b1==atom_index, b1-1==array index
    ai_bonds[b2-1].append(b1)

for i in range(natoms):
    Bgf_Atom.dump_connect(outf, i, ai_bonds[i]) # sort ai_bonds[i] and dump
"""
    if ai_bonds[i]:
        print "%6s%6d" % ( "CONECT", i+1 ),
        #print "{:6}{:6d}".format('CONECT', i+1),
        for ib in ai_bonds[i]:
            if ib == ai_bonds[i][-1]:
                print "%6d" % ib
            else:
                #print "{:5d}".format(ib),
                print "%6d" % ib,
    else:
        print "%6s%6d" % ( "CONECT", i+1 )
"""
outf.write('END\n')

#bonds are not accounted

inf.close()
outf.close()
