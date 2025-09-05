#!/home/joonho/anaconda3/bin/python

import argparse
import re
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from my_data import extract_numeric_block

matplotlib.use('Agg')
def run_calculation(fname, job, icols=[1]):
    with open(fname, 'r') as f:
        lines = f.readlines()
    ### Make 2D
    data2D = extract_numeric_block(lines)
    if len(icols) == 1:
        arr2D = np.array(data2D)[:,icols[0]]
    ### from icols[0] to icols[1]
    elif len(icols) == 2:
        arr2D = np.array(data2D)[:,icols[0]:icols[1]+1]

    print(f"shape {arr2D.shape}")

    if 'fft' in job:
        ### FFT on y-value
        #y = arr2D[:, ic]
        y = arr2D       # only data files in 2D
        #print(f"{y}")
        if job == 'fft':
            fft_y = np.fft.fft(y)
            ### Frequency axis
            n = len(list(y))
            freqs = np.fft.fftfreq(n)
            ### Plot magnitude
            plt.plot(freqs, np.abs(fft_y))
            plt.title("FFT of y-values")
            plt.xlabel("Frequency")
            plt.ylabel("Magnitude")
        elif job == 'ifft':
            fft_y = np.fft.ifftshift(y, axis=0)
            y_real = np.real(fft_y)
            print(y_real)
            np.savetxt("data.dat", y_real, fmt="%.6f", delimiter="\n")
            ### plot
            plt.figure(figsize=(8, 4))
            #plt.plot(y_real[:, 0], label='0.0')
            #plt.plot(y_real[:, 1], label='0.1')
            #plt.plot(y_real[:, 1]-y_real[:, 0], label='diff')
            #plt.title('IFFT Result (Real Part)')
            #plt.xlabel('Index')
            #plt.ylabel('Amplitude')
            #plt.legend()
            #plt.grid(True)
            #plt.tight_layout()
        #plt.show()
        plt.savefig("plot.png")

    else:
        print("Only FFT of 1D is running")
    return 0

def main():
    parser = argparse.ArgumentParser(description='my calculation')
    parser.add_argument('file', help='read data file')
    parser.add_argument('-j', '--job', default='fft', choices=['fft', 'ifft'], help='run fft ifft')
    parser.add_argument('-c', '--col', type=int, nargs='*', help='fft for ith column')
    parser.add_argument('-u','--usage', action='store_true', help="show usage ")
    args = parser.parse_args()
    
    if args.usage:
        print("my_cal.py H_pot.txt -j fft\
            \n\t\tmy_cal.py H_pot.txt -j ifft\
            ")
        sys.exit(0)
    run_calculation(args.file, args.job, args.col) 

if __name__ == '__main__':
    main()
