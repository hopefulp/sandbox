#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import yes_or_no
from mplot2d import mplot_2c

def alldir(d_suff,d_pre,fname,e_limit,Lexe):
    p_dir = os.getcwd()
    lists = os.listdir('.')
    nfile = 0
    l_ene = []
    l_pico= []
    for sdir in lists:
        if sdir.endswith(d_suff):
            nfile += 1
            print(sdir)
            
            #d_id = sdir.split('_')
            #id_sys = d_id[0]
            potfile = p_dir + '/' + sdir + '/' + fname
            i=0
            with open(potfile, "r") as f:
                for line in f:
                    if re.match(" ",line):
                        line.strip()
                        list_l = line.split()
                        if nfile == 1:
                            l_ene.append(float(list_l[1]))
                            l_pico.append(float(list_l[0]))
                        else:
                            l_ene[i] += float(list_l[1])
                        i+=1
                
                        #print(l_ene)


            #cmd = "cat %s " % potfile
            #os.system(cmd)
            #odir = fname[0]+d_rerun
            """ this is for g_energy calculation
            dname = p_dir + '/' + sdir
            os.chdir(dname)
            cmd = "echo 9 | g_energy -f %s.edr -o %s" % (id_sys, "pot")
            if Lexe or yes_or_no("%s ?" % cmd):
                os.system(cmd)
            os.chdir(p_dir)
            """                
            
        #if nfile == 1:
        #    return 0
    print("total number of files: %d, nlen(list): %d" % (nfile, len(l_ene)))
    ave_ene = [ x/float(nfile) for x in l_ene] 
    for x, y in zip(l_pico, ave_ene):
        print("%10.3f%15.5f" % (x, y))

    ### plot S(t)
    numerator = [ x-e_limit for x in ave_ene ]
    den = ave_ene[0] - e_limit
    corr_t = [ x/den for x in numerator ]
    fig_name="CORR" + d_suff
    Lsave=False
    title="CORR_" + d_suff
    xt = "time(ps)"
    yt = "E(kJ/mol)"
    
    mplot_2c(l_pico, corr_t,Lsave,fig_name,title,xt,yt)

    return 0

def main():
    parser = argparse.ArgumentParser(description='execution of all gro file')
    #parser.add_argument('job', choices=["ene"], help='type of input for gromacs')
    parser.add_argument('-s', '--dsuf', help='suffix of directory in search')
    parser.add_argument('-p', '--dpre', help='prefix of searching directory')
    parser.add_argument('-f', '--enefile', help='input energy file')
    parser.add_argument('-e', '--elimit', default=-26.86300, type=float, help='give energy in infinite time')
    # energy in manual is [-26.8630,-363.8034] for chg and full
    parser.add_argument('-r', '--run', action='store_true',help='run without asking')
    args = parser.parse_args()
    
    alldir(args.dsuf,args.dpre,args.enefile,args.elimit,args.run)
    return 0

if __name__=='__main__':
	main()	


