#!/home/jackjack5/epd/bin/python

import argparse
import subprocess
import re

def cmd(id, command, job, freq_option):
    if not id:
        id = subprocess.check_output('whoami')
        idnew = id.rstrip('\n')
    if not command:
        command = 'qstat'
    #print idnew
    str1a = command + ' | grep ' + idnew # + ' ; '
    #str1b = command + ' | grep ' + idnew + ' | grep \' R \' | wc -l '
    #str1 = str1a + str1b
    #print str1
    qstat_sim = subprocess.check_output(str1a, shell=True)

    qstat_list = re.split("\n", qstat_sim)
    print("\n".join(qstat_list))
    file_list=[]
    for line in qstat_list:
        line_list = re.split("\s+", line)
        #print line_list
        if len(line_list) >= 2:
            fname =  line_list[1] + ".out"
            file_list.append(fname)
    #print("\n".join(file_list))
    if job :
        for file in file_list:
            print file,":"
            if job  == "thank":
                subprocess.call(['grep', 'Thank', file ])
            elif job == "freq":
                foption=str(freq_option)
                subprocess.call(['get_freq.pl', file, foption])
            elif job == "efin":
                subprocess.call(['grep', 'Final energy', file])

def main():
    parser = argparse.ArgumentParser(description = 'execute cmd related with qstat')
    parser.add_argument('-p', '--prefix', type = str, help='filename prefix')
    parser.add_argument('-j', '--job', type =str, choices=['thank','freq', 'efin', 'out'])
    parser.add_argument('-f', '--freq_option', type=int, choices=[1, 2], default=1)
    args = parser.parse_args()
    cmd(args.whoami, args.command, args.job, args.freq_option)
    return

if __name__ == '__main__':
    main()
