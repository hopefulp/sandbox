#!/home/jackjack5/epd/bin/python

import argparse
import re
import os

def qsubmit_job(pbs, jobname, ppn, time):
    fout=open("t.csh", 'w')
    with open(pbs, 'r') as f:
        for line in f:          # line has "\n"
            #print line,         # print writes its own "\n"
            if re.search('-N', line):
                strn='#PBS -N ' + jobname + "\n"
                fout.write(strn)
            elif re.search('ppn=', line):
                strn='#PBS -l nodes=1:ppn=' + ppn + "\n"       # type-int cannot be concatenated
                fout.write(strn)
            elif re.search('walltime', line):
                strn='#PBS -l walltime='+time+':00:00'+"\n"
                fout.write(strn)
            else:
                fout.write(line)    # write does not write "\n"
    fout.close()
    str1='qsub t.csh'
    print str1
    x=os.system(str1)
    return x

def main():
    parser = argparse.ArgumentParser(description='change jobname ppn and submit pbs job')
    parser.add_argument('pbsfile', help='basic pbs file')
    parser.add_argument('jobname', help='jobname in pbs file')
    parser.add_argument('ppn', help='number of processor per node in the server')
    parser.add_argument('time', help='run time might be long')
    args = parser.parse_args()
    val = qsubmit_job(args.pbsfile, args.jobname, args.ppn, args.time) 
    if val:
        print "error in qsubmit of job=%s", args.jobname
        exit(3)
    return 

if __name__ == '__main__':
    main()
