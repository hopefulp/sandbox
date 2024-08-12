#!/home/joonho/anaconda3/bin/python

import argparse

coms=[]


def run_amp_sh():
    coms.append("amp_validation.sh -h")
    coms.append("amp_validation.sh Ethylene.extxyz "3 3 3" 0.001 +g | sh")

def run_amp_py():
    coms.append("amp_ene.py Ethylene.extxyz train  -hl 4 4 -el 0.001  -n 5")
    coms.append("amp_ene.py Ethylene.extxyz test  -hl 4 4 -el 0.001  -n 5")
    coms.append("amp_ene.py Ethylene.extxyz md")


def jobs(job, script):
    
    if job == 'ase':
        if script == 'py':
            pass
        elif script == 'sh':
            pass
    elif job == 'amp':
        if script == 'py':
            pass
        elif script == 'sh':
            pass
    else:
        run_amp_sh()
        run_amp_py()
    print(coms)

def main():
    parser = argparse.ArgumentParser(description="display all the works in directory 'py_ai'")
    parser.add_argument('-j', '--job', help="job can be 'ase', 'amp' ")
    parser.add_argument('-s', '--script', help="script can be 'py', 'sh'")
    args = parser.parse_args()

    jobs(args.job, args.script)

if __name__ == "__main__":
    main()
