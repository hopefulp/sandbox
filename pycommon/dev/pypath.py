#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import subprocess
from server_env import home


def run_py(command):
    coml = command.split()

    swhich = f"which {coml.pop(0)}"
    pycommand = subprocess.check_oputput(swhich, shell=True)
    comshell = f'python {pycommand} '
    for com in coml:
        comshell += com + " "
    print(comshell)        
    os.system(comshell)
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="Run python with different HOME in server")
    parser.add_argument('-c','--command', nargs='+',  help="commmand line in shell ")
    args = parser.parse_args()

    run_py(args.command)

if __name__ == "__main__":
    main()
