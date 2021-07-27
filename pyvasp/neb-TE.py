import os, sys, glob

def getdirs():
    dirs = []
    flist = filter(os.path.isdir, glob.glob('*'))
    for fname in flist:
        dirs.append(fname)
    dirs.sort()
    return dirs

def file_read(fname):
    lineinfo = []
    wordinfo = []
    with open(fname) as f:
        for i, l in enumerate(f):
            line = l
            word = line.split()
            lineinfo.append(line)
            wordinfo.append(word)

    return lineinfo, wordinfo

def neb_getenergy(flist):
    FOLDERNAME      = []
    VASPTOTALENERGY = []
    for fname in flist:
        os.chdir(fname)
        if os.path.isfile('OUTCAR'):
            lineinfo, wordinfo = file_read('OUTCAR')
            VASP_E = []
            for i in range(len(lineinfo)):
                if 'y  w' in lineinfo[i]:
                    TE = float(wordinfo[i][-1])
                    VASP_E.append(TE)
                else:
                    pass
            last_E = float(VASP_E[-1])
            FOLDERNAME.append(str(fname))
            VASPTOTALENERGY.append(last_E)
        else:
            pass
        os.chdir('..') 

    for i in range(len(FOLDERNAME)):
        print(FOLDERNAME[i], VASPTOTALENERGY[i])

folder = getdirs()
neb_getenergy(folder)




