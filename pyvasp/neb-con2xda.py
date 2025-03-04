#!/home/joonho/anaconda3/bin/python3

# Min Jong Noh
# Update : 2021 / 03 / 26

'''
This script gathers CONTCAR files in neb calculation (01/02/03 file name...)
and makes XDATCAR file for easy visualization.
'''

import os, sys, glob

usage = 'python3 neb-con2xda.py (image number)'
# ex. python3 neb-con2xda.py 15
n_imag = int(sys.argv[1])
os.system('rm -rf neb-files')

def fileread(fname):
    lineinfo = []
    wordinfo = []
    with open(fname) as f:
        for i, l in enumerate(f):
            line = l
            word = line.split()
            lineinfo.append(line)
            wordinfo.append(word)
    return lineinfo, wordinfo

flist = []
contlist = []
for num in range(n_imag):
    flist.append('%02d' % (num+1))
    contlist.append('CONTCAR_%02d' % (num+1))

os.system('mkdir neb-files')
for i in range(len(flist)):
    os.system('cp %s/CONTCAR neb-files/CONTCAR_%s' % (flist[i], flist[i]))

os.chdir('neb-files')

# first CONTCAR checklist
ori_line, ori_word = fileread('CONTCAR_01')
if ori_word[0] == ori_word[5]:
    VASP_type = 4
else:
    VASP_type = 5

# total number of atoms
natoms = 0
for i in range(len(ori_word[6])):
    natoms += int(ori_word[6][i])

# file writing
XDAT = open('XDATCAR', 'w')

if VASP_type == 4:
    for i in range(7):
        XDAT.write(ori_line[i])
    XDAT.write("Direct configuration=     1\n")
    for i in range(8,8+natoms):
        XDAT.write(ori_line[i])
else: # not working yet
    for i in range(5):
        XDAT.write(ori_line[i])
    XDAT.write(ori_line[0])
    for i in range(7+natoms):
        XDAT.write(ori_line[i])

# CONTCAR file parsing
for j in range(len(contlist)-1):
    con_line, con_word = fileread(contlist[j+1])
    XDAT.write("Direct configuration=     %s \n" % (j+2))
    for i in range(8, 8+natoms):
        XDAT.write(con_line[i])
XDAT.close()
