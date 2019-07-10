#!/home/joonho/anaconda3/bin/python

import os
import argparse

d_search=['/home/joonho/bin', '/home/joonho/sandbox_gl']
d_backup='/shared/share_win/backup_script'

def search_dir(d,Lrun):
    ### change directory to current directory
    os.chdir(d)
    print(f"####.... enter {d}")
    files = os.listdir('.')
    for f in files:
        if os.path.isdir(f):
            ### enter into sub-directory
            os.chdir(f)
            cwd = os.getcwd()
            search_dir(cwd, Lrun)
            ### should come back to the present directory after finish search sub-directory
            os.chdir('..')
        else:
            fname = f"{d}/{f}"
            if os.access(fname, os.X_OK):
                if os.path.islink(fname):
                    com = f"rsync -avz -L {fname} {d_backup}"
                else:
                    com = f"rsync -avz {fname} {d_backup}"
                if Lrun:
                    os.system(com)
                else:
                    print(com)
    print(f"####.... exit {d}")
def jobs(dir_input, Lrun):
    ### gather directories to be searched
    if dir_input:
        for d in dir_input:
            d_search.append(d)

    for d in d_search:
        if os.path.isdir(d) and os.path.isdir(d_backup):
            search_dir(d, Lrun)
        else:
            print(f"directory deos not exist: {d} or {d_backup}")
            exit(0)
    if not Lrun:
        print("Use -r to run")
    return 0

def main():
    parser = argparse.ArgumentParser(description="backup all the scripts  ")
    parser.add_argument('-d','--directory', nargs='*', help="add more directories which contains scripts ")
    #parser.add_argument('-j','--job', help="present class details ")
    #parser.add_argument('-js','--specify', choices=['qcmo','nbo','eda'], help="present class details ")
    parser.add_argument('-r','--run', action='store_true', help="execute command")
    args = parser.parse_args()

    jobs(args.directory, args.run)

if __name__ == "__main__":
    main()
