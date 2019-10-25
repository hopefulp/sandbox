#!/usr/bin/python

import os
import argparse

def allvasp():
    lists = os.listdir('.')
    for dir in lists:
        if not os.path.isdir(dir):
            continue
        os.chdir(dir)
        print('now on' + dir)
        if os.path.isfile('POSCAR') and os.path.isfile('POTCAR') and not os.path.isfile('OUTCAR'):
            os.chdir(dir)
            os.system('echo "here is", dir')
        os.chdir('..')

def main():
    parser = argparse.ArgumentParser(description='execution of all vaspallallall')
    args=parser.parse_args()
    allvasp()
    return


if __name__ == '__main__':
    main()
