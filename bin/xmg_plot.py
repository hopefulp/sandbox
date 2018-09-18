#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys


def main():
    """ This program is used for 
    1. all electron
        input: xmg_plot.py 1st/si-lda-ae -j ae
        xmgrace  -nxy 1st/si-lda-ae.wfc -p $gracehome/Siae3.par &
    2. ps - now default
        xmg_plot.py 1st/si-lda
        xmgrace  -nxy 1st/si-ldaps.wfc -p $gracehome/Si_pp.par &
    3. ld (ld or ld2-including ps)
        ld
            xmg_plot.py si-LDcheck -j ld
            xmgrace  -nxy si-LDcheck.dlog -p $gracehome/Si_ld.par &
        ld2            
            xmg_plot.py si-LDcheck -j ld2
            xmgrace  -nxy si-LDcheck.dlog -nxy si-LDcheckps.dlog -p $gracehome/Si_ld.par &
        ld2 -c (2 directory)
            xmg_plot.py si-TS -n 5 -j ld2 -c 2nd
            xmgrace  -nxy si-TS5.dlog -nxy si-TS5ps.dlog -nxy 2nd/si-TS5ps.dlog -p $gracehome/Si_ld2f.par &
    4. p for ionic state
        xmg_plot.py si-TS -j val -b p -n 1 2 3 4 5
        xmgrace  -block si-TS1ps.wfc -bxy 1:3 -block si-TS2ps.wfc -bxy 1:3 -block si-TS3ps.wfc -bxy 1:3 -block si-TS4ps.wfc -bxy 1:3 -block si-TS5ps.wfc -bxy 1:3 -p $gracehome/Si_plus.par &
        for s
        xmg_plot.py si-TS -j val -b s -n 1 2 3 4 5 -p Si_plus1.par
        
    """
    parser = argparse.ArgumentParser(description="xmgrace running")
    parser.add_argument('fname', help='give filename prefix')
    parser.add_argument('-j', '--job', default='ps', choices=['ps', 'ae', 'ld', 'ld2', 'val'], help='graph kinds of pseudo potential or logarithmic derivative')
    #parser.add_argument(
    parser.add_argument('-n', '--nfile', default=[''], nargs='*', help='number for serial ionic state where charge is n-1')
    parser.add_argument('-d', '--direct', help='two files are drawn in one figure for different directory')
    parser.add_argument('-f', '--fname2', help='add file for comparison')
    parser.add_argument('-b', '--block', choices=['s', 'p', 'd'], help='draw selected column for several files')
    parser.add_argument('-p', '--par', help='input *.par file') 

    args = parser.parse_args()
    nlist= args.nfile
    files= []
    files2=[]
    if args.job == 'ae':
        suff = ".wfc"
        pfile = "$gracehome/Siae3.par"
    elif args.job == 'ps' or args.job == 'val':
        suff = "ps.wfc"
        if args.job == 'ps':
            pfile = "$gracehome/Si_pp.par"
        elif args.job == 'val':
            pfile = "$gracehome/Si_plus.par"
    ### for LD, draw ae and ps        
    elif args.job == 'ld' or args.job == 'ld2':
        suff = ".dlog"
        pfile = "$gracehome/Si_ld5.par"
        if args.direct:
            pfile = "$gracehome/Si_ld2f.par"
        if args.job == 'ld2':
            suff2 = "ps.dlog"
    else:
        print("error in job\n")
        sys.exit(1)

    if args.par:
        pfile = "$gracehome/" + args.par

    for n in args.nfile:
        fname = args.fname + str(n) + suff
        files.append(fname)
        if args.job == 'ld2':
            fname2 = args.fname + str(n) + suff2
            files2.append(fname2)
    print(files)

    cmd1 = 'xmgrace '
    
    if args.block: # for p or s separately in several files
        cmd2 = ''
        if args.block == 'p':
            ae_col = 3
        elif args.block == 's':
            ae_col = 2
        else:
            print("d-orbital is not included in valence yet")
            sys.exit(10)
        for f in files:
            cmd2 += ' -block ' + f + ' -bxy 1:' + str(ae_col)
        if pfile:
            cmd3 = ' -p ' + pfile
        cmd = cmd1 + cmd2 + cmd3 + ' &'        
        print(cmd)

        os.system(cmd)
    else:
        # in case two files are drawn (ld2)
        if len(files) == len(files2):
            for (f, f2) in zip(files, files2):
                cmd2 = ' -nxy ' + f
                if args.job == 'ld2': # for just 1 file
                    cmd2 += ' -nxy ' + f2
                if args.direct:
                    fname3 = args.direct + '/' + f2
                    cmd2 += ' -nxy ' + fname3

                if pfile:
                    cmd3 = ' -p ' + pfile
                cmd = cmd1 + cmd2 + cmd3 + '  &'        
                print(cmd)
                os.system(cmd)
        else:               
            for f in files:
                cmd2 = ' -nxy ' + f
                if args.job == 'ld2': # for just 1 file
                    cmd2 += ' -nxy ' + fname2
                if args.direct:
                    fname3 = args.direct + '/' + fname2
                    cmd2 += ' -nxy ' + fname3

                if pfile:
                    cmd3 = ' -p ' + pfile
                cmd = cmd1 + cmd2 + cmd3 + '  &'        
                print(cmd)
                os.system(cmd)
    return 0            

if __name__ == "__main__":
    main()
