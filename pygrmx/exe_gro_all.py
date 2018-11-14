#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
from common import yes_or_no

def allgro(ftype,d_suff,ftop,fmdp,fin,d_rerun,o_trj,Lexe):
    home = os.getenv('HOME')
    p_dir = os.getcwd()
    lists = os.listdir('.')
    pbsfile=home + '/sandbox_gl/pypbs/sge_mdrun.tcsh'
    nfile = 0
    for fgro in lists:
        if fgro.endswith(ftype):
            nfile += 1
            fname = fgro.split('.')
            ndir = fname[0]+d_suff
            odir = fname[0]+d_rerun
            new_dir = p_dir + '/' + ndir
            old_dir = p_dir + '/' + odir
            os.mkdir(new_dir)
            cmd = 'cp %s %s %s ' % (fgro, ftop, fmdp)
            if fin:
                cmd += ' %s ' % fin
            cmd += ' %s ' % new_dir
            print(cmd)
            if Lexe or yes_or_no('%s ?' % cmd):
                os.system(cmd)
            if d_rerun:
                f_rerun = odir + '/' + fname[0] + '.trr'
                target = ndir + '/excited.trr'
                cmd = 'cp %s %s' % ( f_rerun, target )
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
            os.chdir(p_dir)
            cmd = 'mv %s tmp' % fgro
            os.system(cmd)

        #if nfile == 1:
        #    return 0
    return 0

def main():
    parser = argparse.ArgumentParser(description='execution of all gro file')
    parser.add_argument('ftype', default='gro', help='type of input for gromacs')
    parser.add_argument('dsuff', help='add suffix of new directory')
    parser.add_argument('-t', '--topology', help='input topology file')
    parser.add_argument('-p', '--mdparam', help='md parameter file')
    parser.add_argument('-i', '--include', help='include itp file file')
    parser.add_argument('-r', '--rerund', help='include trr file for rerun')
    parser.add_argument('-j', '--ftrj', default='excited.trr',help='given trajectory file for rerun')
    parser.add_argument('-e', '--exe', action='store_true',help='run without asking')
    args = parser.parse_args()
    
    allgro(args.ftype,args.dsuff,args.topology,args.mdparam,args.include,args.rerund,args.ftrj,args.exe)
    return 0

if __name__=='__main__':
	main()	


