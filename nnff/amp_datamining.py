#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys
import numpy as np

from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp.descriptor.analysis import FingerprintPlot
from amp.utilities import hash_images, get_hash
from amp import Amp

from ase import io
from myplot2D import barplot_y
import amp_util

def data_fp(inf, job, sub_command,  i_image, iatom, Lgraph, xtitle, ytitle, title):
    ### 
    str_index = amp_util.list2str(i_image)
    images = io.read(inf, index=str_index)
    nimages = len(images)
    natoms  = len(images[0])
    ### images changes into hash dictionary
    ### fp needs connection to Amp; but Gs is random
    ### calculate_fingerprints makes descriptor has fingerprints, format: 'O', [fp lists]
    ### fps: natoms, fps[0]: 1st atom tuple of ('O', feature-vec), fps[0][1]: feature-vector
    hash_dict = hash_images(images)
    try:
        calc = amp_util.get_amppot()
    except FileNotFoundError:
        calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(4))
    calc.descriptor.calculate_fingerprints(hash_dict)
    ### to print here
    hash0 = get_hash(images[0])
    fps = calc.descriptor.fingerprints[hash0]
    print(f"nimages {nimages}, natoms {natoms} nfeatures(fp) {len(fps[0][1])}")
    for image in hash_dict:
        fps = calc.descriptor.fingerprints[image]
        if not iatom == None:
            fp = fps[iatom]    
        #for fp in fps:
            print(fp[0], " ".join("%7.4f" % x for x in fp[1]))
            
    if Lgraph:
        fpplot = FingerprintPlot(calc)
        fpplot._calc.descriptor.calculate_fingerprints(hash_dict)
        print(f"{len(i_image)}: {iatom}")
        if len(i_image) == 1 and iatom != None:
            fingerprints = fpplot._calc.descriptor.fingerprints[list(hash_dict.keys())[0]]
            ts = f"{fingerprints[iatom][0]}-{iatom}/{i_image}-image"
            if title:
                ts += f": {title}"
            ### [0] for atom species [1] for fp
            barplot_y(fingerprints[iatom][1], xlabel=xtitle, ylabel=ytitle, title=ts)
            print(fingerprints[iatom][1])
            print("print an atom")
        else:
            fpplot(hash_dict)   # FingerprintPlot has __call__(self, images)
            print(f"The fingerprint range of {nimages} images from {i_image[0]} were presented in 'fingerprints.pdf'")
    return

def main():
    parser = argparse.ArgumentParser(description='Data analysis ', prefix_chars='-+/')
    parser.add_argument('-f', '--infile', default='OUTCAR', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    ### tef for test w. force
    parser.add_argument('-j', '--job', default='fp', choices=['fp'], help='job option:fp-fingerprint')
    parser.add_argument('-s', '--sub_command', default='atom', choices=['atom'], help='subcommand for atom')
    #parser.add_argument('-p', '--pot', help="input amp potential")
    ### DATA group
    data_group = parser.add_argument_group(title='Data')
    data_group.add_argument('-im', '--iimage', nargs='+', default=[0], type=int, help='index of image')
    data_group.add_argument('-ia', '--iatom', type=int, help='index of atom')
    ### PLOT group
    plot_group = parser.add_argument_group(title='Plot')
    plot_group.add_argument('+g', action="store_true", help='if val default is False, otherwise True')
    plot_group.add_argument('-xt', '--xtitle', default='fingerprint', help='xtitle in plot')
    plot_group.add_argument('-yt', '--ytitle', default='f(fp)', help='ytitle in plot')
    plot_group.add_argument('-t', '--title', help='ytitle in plot')

    args = parser.parse_args()

    pwd = os.getcwd()
    if args.infile:
        fname=args.infile
    else:
        fname=amp_util.get_inputfile(pwd)
    if not fname:
        sys.exit(1)

    data_fp(args.infile, args.job, args.sub_command, args.iimage, args.iatom, args.g, args.xtitle, args.ytitle, args.title)
    return

if __name__ == '__main__':
    main()

