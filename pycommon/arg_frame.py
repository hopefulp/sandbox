#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import yes_or_no
from mplot2d import mplot_2c

def f_sub(f1,f2,Lexe):
    p_dir = os.getcwd()
    lists   = os.listdir('.')
    nfile   = 0
    l_ene   = []
    l_pico  = []
    l_pico2 = []
    for sdir in lists:
        if sdir.endswith(d_suff):
            nfile += 1
            
            d_id = sdir.split('_')
            sdir2 = d_id[0]+'_'+d_suff2
            potfile = p_dir + '/' + sdir + '/' + fname
            potfile2= p_dir + '/' + sdir2+ '/' + fname
            print("%s %s" % (sdir, sdir2))
            i=0
            with open(potfile, "r") as f, open(potfile2, "r") as g:
                for line, line2 in zip(f, g):
                    if re.match(" ",line):
                        line.strip()
                        line2.strip()
                        list_l = line.split()
                        list_l2= line2.split()
                        if nfile == 1:
                            #l_ene.append(float(list_l[1])) # for single file E
                            del_ene = float(list_l[1]) - float(list_l2[1])
                            l_ene.append(del_ene)
                            l_pico.append(float(list_l[0]))
                        else:
                            #l_ene[i] += float(list_l[1])   # for single file E
                            del_ene = float(list_l[1]) - float(list_l2[1])
                            l_ene[i]+= del_ene
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
    outf = fig_name + '.dat'
    Lsave=False
    title="CORR_" + d_suff
    xt = "time(ps)"
    yt = "E(kJ/mol)"
    
    mplot_2c(l_pico, corr_t,Lsave,fig_name,title,xt,yt)
    print("outfile: %s" % outf)
    with open(outf, "w") as f:
        for x, y in zip(l_pico, corr_t):
            print("%10.3f%15.5f" % (x, y))
            f.write("%10.3f%15.5f" % (x, y))


    return 0

def main():
    parser = argparse.ArgumentParser(description='execution of all gro file')
    #parser.add_argument('job', choices=["ene"], help='type of input for gromacs')
    parser.add_argument('-f1', '--file1', help='suffix of directory in search')
    parser.add_argument('-f2', '--file2', help='suffix of directory in search')
    parser.add_argument('-r', '--run', action='store_true',help='run without asking')
    args = parser.parse_args()
    
    f_sub(args.file1,args.file2,args.run)
    return 0

if __name__=='__main__':
	main()	


