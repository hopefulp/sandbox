#!/home/jackjack5/epd/bin/python
### 2017.06.22 dev1

import argparse
import os
import re
import sys

bin='/qcfs/joonho/bin/'
#cmd="qgeo_mo.pl"
PPN={'kdft':12, 'psi':16, 'rho': 16, 'mu':32}

def all_qchem(remf, pbsfile):
    nmol=0
    ### Obtain ppn
    host = os.getenv('HOST')
    if host == "" or host not in PPN.keys():
        print "host is not identified"
        sys.exit(10)
    else:
        ppn=PPN[host]
    ### Scan pwd to find .mol file
    list = os.listdir('.')
    for fname in list:
        # if not .mol file, skip
        if not re.search('\.mol$', fname):  # if not \. in filename, skip
            #print "not molfile: ", fname
            continue
        #else:           # re practice
        #    print " yes : ", fname
        #    continue
        lfname=re.split("\.",fname)
        # if there is not only one \., skip
        if lfname[1] != "mol":
            print "Error in reading .mol filename in ", fname
            continue
        outf=lfname[0]+'.in'
        qchem_outf=lfname[0]+'.out'
        # if there is already qchem outfile, skip
        if os.path.isfile(qchem_outf):
            continue
        nmol += 1
        # get mol from qchem.out
            #str1 = cmd + " " + fname 
            #if re.search('ts', fname):
            #    str1 += " 3"
            #print str1
            #os.system(str1)
        ### Make qchem input file
        #if re.search('spts', fname):
        #    suff='spts'
        if re.search('ts', fname):
            suff='ts'
        elif re.search('irc', fname):
            suff='irc'
        else:
            suff='opt'
        remfile = remf + suff
        # if rem. file does not exist, exit w. error message
        if not os.path.isfile(remfile):
            print remfile, " does not exist"
            sys.exit(3)
        str1 = 'cat ' + fname +' '+ remfile + ' > ' + outf
        print str1
        os.system(str1)
        ### Submit job
        str2 = 'qsub_qc_changeline.py ' + pbsfile + ' '+ lfname[0] + ' ' + str(ppn)
        print str2
        os.system(str2)
        # cannot move mol file before job in queue is running
        #if nmol == 1:
        #    if not os.path.isdir('Jobdone'):
        #        os.system('mkdir Jobdone')
        #os.system('mv fname Jobdone')

    if nmol != 0:
        print "%d job was submitted" % nmol
    else:
        print "no job was submitted"
        sys.exit(6)
    return


def main():
    parser = argparse.ArgumentParser(description='execute qchem for all .mol files if not .outfile')
    parser.add_argument('rem', type=str, help='rem file for qchem input')
    parser.add_argument('pbs', type=str, help='pbs file for qchem submit')
    args = parser.parse_args()
    all_qchem(args.rem, args.pbs)
    return    

if __name__ == '__main__':
    main()


