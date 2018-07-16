#!/usr/bin/python

import argparse
import os
import Env_msg      # Machine -> Env_msg
import re

def vasp_jobs(dirs, run):

    host = Env_msg._HOST
    pbs = Env_msg.host_machine.vasp_pbs
    pbs_name = Env_msg.pbs_vname
    commands=[]
    if not dirs:
        print "input directory to run out of directory"
        exit(0)
    elif re.match("all", dirs[0], re.IGNORECASE):
        lists = os.listdir('.')
        for dir in lists:
            if not os.path.isdir(dir):
                continue
            os.chdir(dir)
            print 'now on '+dir
            if os.path.isfile('INCAR') and os.path.isfile('POSCAR') and not os.path.isfile('OUTCAR'):
                os.system('csh /qcfs/jackjack5/vasp/vaspenv.sh')
            os.chdir('..')
    else:
        for dir in dirs:
            if not os.path.isdir(dir):
                print "%f is not directory" % dir
                exit(11)
            print "copy pbs from ..: ", pbs
            cmd = 'cp ' + pbs + ' .'
            os.system(cmd)
            cmd = "sed -i 's/#PBS -N/#PBS -N " + dir + "/' " + pbs_name
            print cmd
            os.system(cmd)
            cmd = 'qsub ' + pbs_name
            print cmd
            if run:
                os.system(cmd)
        if not run:
            print "to run use option -r"
            exit(1)

    return 0

def main():
    parser = argparse.ArgumentParser(description='execution of vasp for directories')
    parser.add_argument('dir', nargs='*', help='input directories or all for all the directories')
    parser.add_argument('-r', '--run', action='store_true', help='run qsub')
    args = parser.parse_args()
    vasp_jobs(args.dir, args.run)
    
    return 0

if __name__=='__main__':
	main()	


