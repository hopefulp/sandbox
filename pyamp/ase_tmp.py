#!/home/joonho/anaconda3/bin/python

import argparse

import ase.io

import numpy as np

def anal_traj(traj_file, job):

    if job == "traj":

        traj = ase.io.Trajectory(traj_file, "r")
        


def main():
    parser = argparse.ArgumentParser(description='analize ASE trajectory')
    parser.add_argument('fin', help='ASE traj file')
    parser.add_argument('job', default='traj', choices=['traj'], help='job option: "traj"')
    args = parser.parse_args()

    anal_traj(args.fin, args.job)
    return

if __name__ == '__main__':
    main()

