#!/usr/bin/python
### 2017.06.22 dev1
### 2018.04.23 dev2

import argparse
import os
import re
import sys
from Machine import *       # PPNs, PBS_QChem_scripts
from common import *

rem_files=[]
all_files=[]
pbsfile_dir = os.getenv('MODULE_SOURCE') + '/pypbs/'

### get remfile for each mol file
def get_rem_files():
    global all_files

    if not all_files:
        all_files = os.listdir('.')
    if not rem_files:
        for fname in all_files:
            if re.match('rem', fname):
                rem_files.append(fname)

def run_qchem_all(m_tag, patterns, Rem_file, isrun, file_type, pbs_long, ncore):
    global all_files
    nmol=0
    ### Obtain ppn and pbsfile
    host = os.getenv('HOST')
    pwd = os.getcwd()
    if host == "" or host not in PPNs.keys():
        print "host is not identified"
        sys.exit(10)
    else:
        if ncore:
            ppn = ncore
        else:
            ppn=PPNs[host]
        pbsfile = pbsfile_dir + PBS_QChem_scripts[host]
    print "PBS file :: ", pbsfile 
    ### Scan pwd to get file list to be run [mol|in|specified]
    ## if File_list is not specified
    l_file=[]
    # in case file list
    if m_tag == 'i':
        l_file = patterns[:]
    # otherwise get files by patterns
    else:
        l_file = get_files_pattern(m_tag, patterns, pwd)
    for file in l_file:
        lfile_name=re.split("\.",file)
        if len(lfile_name) >= 3:
            print "Filename with multiple period:: Error "
            exit(10)
        qchem_in=lfile_name[0]+'.in'
        qchem_outf=lfile_name[0]+'.out'
        if lfile_name[1] == 'in':
            file_type = 'in'
        ### SKIP if there is out file
        if os.path.isfile(qchem_outf):
            #if File_list:
            print("As for %s, there exists .out file" % file)
            continue
        else:
            nmol += 1

        ### LOOKing for rem file
        if file_type == "mol" and nmol == 1 :
            if Rem_file:
                rem_file = Rem_file
            else:
                rem_prefix = ["rem\."]
                rem_files = get_files_prefix(rem_prefix, pwd)
                if len(rem_files) == 1:
                    rem_file = rem_files[0]
                else:
                    #print rem_files
                    if re.search('irc', file, re.IGNORECASE):
                        Rfiles = [ s for s in rem_files if 'irc' in s]   # [   ] returns list
                    elif re.search('ts', file, re.IGNORECASE):
                        Rfiles = [ s for s in rem_files if 'spts' in s]
                    else:
                        Rfiles = [ s for s in rem_files if 'opt' in s]
                    rem_file = Rfiles[0]
                    q = "will you use %s file" % rem_file
                    if not yes_or_no(q):
                        print "specify rem file among -- ", Rfiles
                        exit(20)
        ### RUN if it's ready
        ## Make in file in case of mol file
        if file_type == "mol":
            str1 = 'cat ' + file +' '+ rem_file + ' > ' + qchem_in
            print str1
            if isrun :
                os.system(str1)
        ## Submit job : pbs file read jobname.in
        str2 = 'qsub_changeline_qc.py ' + pbsfile + ' '+ lfile_name[0] + ' ' + str(ppn) + ' ' + pbs_long
        print str2
        if isrun:
            os.system(str2)
            
    if nmol != 0:
        if isrun:
            print "%d job was submitted" % nmol
        else:
            print "%d job will be submitted" % nmol
            print "to run use '-e' option"
    else:
        print "no job was submitted"
        sys.exit(6)
    return


def main():
    parser = argparse.ArgumentParser(description='execute qchem for all .mol files if not .out or .in')
    #parser.add_argument('-p', '--pbs', type=str,  help='pbs file for qchem submit')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-i', '--infile',  nargs='+', help='run input files in group 2')
    group.add_argument('-p', '--prefix', nargs='+', help="prefix exclusive to suffix and search")
    group.add_argument('-s', '--suffix', nargs='+', help="suffix exclusive to prefix and search")
    group.add_argument('-m', '--matching', nargs='+', help="matching(re.search) exclusive to prefix and suffix")
    parser.add_argument('-t', '--type', type=str, default='mol', choices=['mol', 'in'])
    parser.add_argument('-e', '--exe', action='store_true', help='run job or just display')
    parser.add_argument('-r', '--rem', help='rem file : rem.file for full name')
    parser.add_argument('-l', '--pbs_long', default='120', choices=['240', '480'], help='run pbs in long time')
    parser.add_argument('-c', '--core', type=int, help='number of core process')
    args = parser.parse_args()

    if args.infile:
        matching=args.infile
        m_tag = 'i'
    elif args.prefix:
        matching=args.prefix
        m_tag = 'p'
    elif args.suffix:
        matching=args.suffix
        m_tag = 's'
    elif args.matching:
        matching=args.matching
        m_tag = 'm'
    else:
        print "matching should be given"
        return 1
    
    run_qchem_all(m_tag, matching, args.rem, args.exe, args.type, args.pbs_long, args.core)
    return    

if __name__ == '__main__':
    main()


