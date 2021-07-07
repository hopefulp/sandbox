#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys
import numpy as np

from amp import Amp
from amp.descriptor.analysis import FingerprintPlot
from amp.utilities import hash_images, get_hash

from ase import io
from myplot2D import barplot_y

def amp_analysis(inf, pot, job, sub_command,  i_image, iatom, xtitle, ytitle, title):
    ### parameter file, basically it doesn't need calc if fp was already there
    calc = Amp.load(pot)
    fpplot = FingerprintPlot(calc)
    if len(i_image) == 1:
        sindex = "%i:%i" % (i_image[0], i_image[0]+1)                 # get only one image
        nimage = 1
    else:
        sindex = "%i:%i" % (i_image[0], i_image[1])
        nimage = i_image[1] - i_image[0]
    images = io.read(inf, index=sindex)
    images = hash_images(images)
    fpplot._calc.descriptor.calculate_fingerprints(images)  # add fingerprints in fpplot._calc.descriptor
    if len(i_image) == 1 and iatom:
        fingerprints = fpplot._calc.descriptor.fingerprints[list(images.keys())[0]]
        ts = f"{fingerprints[iatom][0]}-{iatom}/{i_image}-image"
        if title:
            ts += f": {title}"
        barplot_y(fingerprints[iatom][1], xlabel=xtitle, ylabel=ytitle, title=ts)
    else:
        fpplot(images)
        print(f"The fingerprint range of {nimage} images from {i_image[0]} were presented in 'fingerprints.pdf'")
    return

def find_inputfile(pwd, job):
    if os.path.isfile('OUTCAR'):
        fname = 'OUTCAR'
    else:
        for f in os.listdir(pwd):
            if f.endswith('extxyz'):
                fname = f
                break
    quest = f"file {fname} will be read with job='{job}' ?"
    if yes_or_no(quest):
        return fname
    else:
        return 0
    

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz, OUTCAR: validation is removed ', prefix_chars='-+/')
    parser.add_argument('-f', '--infile', default='OUTCAR', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    ### tef for test w. force
    parser.add_argument('-j', '--job', default='fp', choices=['fp'], help='job option:fp-fingerprint')
    parser.add_argument('-s', '--sub_command', default='atom', choices=['atom'], help='subcommand for atom')
    parser.add_argument('-p', '--pot', help="input amp potential")
    ### DATA group
    data_group = parser.add_argument_group(title='Data')
    data_group.add_argument('-im', '--iimage', nargs='+', type=int, help='index of image')
    data_group.add_argument('-ia', '--iatom', type=int, help='index of atom')
    ### PLOT group
    plot_group = parser.add_argument_group(title='Plot')
    plot_group.add_argument('-xt', '--xtitle', default='fingerprint', help='xtitle in plot')
    plot_group.add_argument('-yt', '--ytitle', default='f(fp)', help='ytitle in plot')
    plot_group.add_argument('-t', '--title', help='ytitle in plot')

    args = parser.parse_args()

    pwd = os.getcwd()
    if args.infile:
        fname=args.infile
    else:
        fname=find_inputfile(pwd, args.job)
    if not fname:
        sys.exit(1)
    amp_analysis(args.infile, args.pot, args.job, args.sub_command, args.iimage, args.iatom, args.xtitle, args.ytitle, args.title)
    return

if __name__ == '__main__':
    main()

