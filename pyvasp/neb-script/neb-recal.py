import os, sys, glob

def getdirs():
    dirs = []
    flist = filter(os.path.isdir, glob.glob('*'))
    for fname in flist:
        dirs.append(fname)
    dirs.sort()
    return dirs

def nebnext(flist):
    for fname in flist:
        os.chdir(fname)
        os.system('clearvasp.sh')
        os.system('mv CONTCAR POSCAR')
        os.chdir('..')
    print('The NEB files are successfully transfered')

folder = getdirs()
nebnext(folder)




