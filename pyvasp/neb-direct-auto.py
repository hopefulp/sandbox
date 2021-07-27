import os, sys, glob

def getdirs():
    dirs = []
    flist = filter(os.path.isdir, glob.glob('*'))
    for fname in flist:
        dirs.append(fname)
    dirs.sort()
    return dirs

def posdirect(flist):
    for fname in flist:
        os.chdir(fname)
        os.system('python3 ~/bin/neb-direct.py')
        os.chdir('..')
    print('The POSCAR files are successfully changed to direct coordinates')

folder = getdirs()
posdirect(folder)




