#!/home/joonho/anaconda3/bin/python
import time
import os
import argparse


def command(lcom,stime):

    comm = " ".join(lcom) 
    print(comm)

    while True:
        os.system(comm)
        time.sleep(stime)

def main():
    parser = argparse.ArgumentParser(description = 'run command continuously')
    parser.add_argument('command', nargs='+',  help='get command line as list')
    parser.add_argument('-t', '--stime', default=10, type=int, help='sleep time for module time.sleep')
    args = parser.parse_args()
    command(args.command, args.stime)

if __name__ == '__main__':
    main()

    
