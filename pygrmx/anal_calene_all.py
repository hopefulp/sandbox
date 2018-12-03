#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
from common import yes_or_no

def alldir(ftype,d_suff,d_pre,Lexe):
    p_dir = os.getcwd()
    lists = os.listdir('.')
    nfile = 0
    for sdir in lists:
        if sdir.startswith(d_pre):
            nfile += 1
            print(sdir)
            
            d_id = sdir.split('_')
            id_sys = d_id[0]

            #odir = fname[0]+d_rerun

            dname = p_dir + '/' + sdir
            os.chdir(dname)
            cmd = "echo 9 | g_energy -f %s.edr -o %s" % (id_sys, "pot")
            if Lexe or yes_or_no("%s ?" % cmd):
                os.system(cmd)
            """
            os.mkdir(new_dir)
            cmd = 'cp %s %s %s ' % (fgro, ftop, fmdp)
            if fin:
                cmd += ' %s ' % fin
            cmd += ' %s ' % new_dir
            print(cmd)
            if Lexe or yes_or_no('%s ?' % cmd):
                os.system(cmd)
            os.chdir(new_dir)
            cmd = 'grompp -o %s.tpr -f %s -c %s -p %s' % (fname[0],fmdp,fgro,ftop)
            #print(cmd)
            if Lexe or yes_or_no('%s ?' % cmd):
                os.system(cmd)
            #cmd = 'qsub -N %s -v tpr=%s %s' % (ndir,fname[0],pbsfile)
            if d_rerun:
                job = "rerun"
                cmd = 'qsub -N %s -v tpr=%s -v job=%s -v trj=%s %s' % (ndir,fname[0],job,ftrj,pbsfile)
            else:
                job = "md"
                cmd = 'qsub -N %s -v tpr=%s -v job=%s  %s' % (ndir,fname[0],job,pbsfile)
            if Lexe or yes_or_no('%s ?' % cmd):
                os.system(cmd)
            cmd = 'mv %s tmp' % fgro
            os.system(cmd)
            """
            os.chdir(p_dir)
            
        #if nfile == 1:
        #    return 0
    return 0

def main():
    parser = argparse.ArgumentParser(description='execution of all gro file')
    parser.add_argument('job', choices=["ene"], help='type of input for gromacs')
    parser.add_argument('-s', '--dsuf', help='suffix of directory in search')
    parser.add_argument('-p', '--dpre', help='prefix of searching directory')
    parser.add_argument('-t', '--topology', help='input topology file')
    parser.add_argument('-i', '--include', help='include itp file file')
    parser.add_argument('-r', '--rerund', help='include trr file for rerun')
    parser.add_argument('-j', '--ftrj', default='excited.trr',help='given trajectory file for rerun')
    parser.add_argument('-e', '--exe', action='store_true',help='run without asking')
    args = parser.parse_args()
    
    alldir(args.job,args.dsuf,args.dpre,args.exe)
    return 0

if __name__=='__main__':
	main()	


