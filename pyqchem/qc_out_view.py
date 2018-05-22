#!/home/jackjack5/epd/bin/python

import argparse
import subprocess
import re
import os

File_object1 = 'recent files'
File_object2 = 'selected files'

def cmd(obj,obj_arg, job, freq_option):
    file_list=[]
    if obj == File_object1:
        num = int(obj_arg)
        line = subprocess.Popen('ls -dlt *.out', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = line.communicate()
        #print out
        lines = re.split("\n", out) 
        for i in range(num):
            line_ele = re.split("\s+", lines[i])
            fname = line_ele[-1]
            #print fname
            file_list.append(fname)
    elif obj == File_object2:
        f_pre = obj_arg
        files = os.listdir('.')
        for fname in files:
            if re.match(f_pre, fname) and re.search('out$', fname):
                #print fname
                file_list.append(fname)
    if job :
        for file in file_list:
            print file,":"
            if job  == "thank":
                subprocess.call(['grep', 'Thank', file ])
            elif job == "freq":
                foption=str(freq_option)
                subprocess.call(['get_freq.pl', file, foption])
            elif job == "efin" or job == "gfin":
                subprocess.call(['grep', 'Final energy', file])

def main():
    parser = argparse.ArgumentParser(description = 'check recent or specified Q-Chem outfile')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-n', '--number', type=int, help='number of Q-Chem out-files')
    group.add_argument('-p', '--prefix', type=str, help='prefix of Q-Chem out-files')
    parser.add_argument('job', type =str, choices=['thank','freq', 'efin','gfin', 'out'])
    parser.add_argument('-f', '--freq_option', type=int, choices=[1, 2], default=1)
    args = parser.parse_args()
    # treat exclusive options
    if args.prefix:
        object = File_object2
        obj_arg = args.prefix
    else :
        object = File_object1
        if args.number:
            obj_arg = args.number
        else:
            obj_arg = 5         # default for number of recent files
    cmd(object, obj_arg, args.job, args.freq_option)
    return

if __name__ == '__main__':
    main()
