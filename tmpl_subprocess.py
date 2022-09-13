#!/home/joonho/anaconda3/bin/python

import argparse
import subprocess
import re

def shell_command(bash, partition):
    popen = subprocess.Popen(bash, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdoutdata, sterrdat = popen.communicate()

    data = stdoutdata.decode('utf-8')

    print(data)

    return 1

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('shell', default='pestat', choices=['pestat'], help='get input from shell command')
    parser.add_argument('-part', '--partition', default='all', help='partition number')
    args = parser.parse_args()
    
    shell_command(args.shell, args.partition) 

if __name__ == '__main__':
    main()
