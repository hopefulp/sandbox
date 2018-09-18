#!/usr/bin/python

import numpy as np
import random
import argparse
import sys
import operator as op
from My_arith import *

def make_spin(L, s_conf):
    M2d=[]
    if s_conf == 'random':
        up_prob = 0.5
    elif s_conf == 'up':
        up_prob = 1.0
    elif s_conf == 'down':
        up_prob = 0.0

    for i in range(L):
        row_vec=[]              # row vector in M2d
        for j in range(L):
            if up_prob == 1.0:
                row_vec.append(1)
            elif up_prob == 0.0:
                row_vec.append(-1)
            else:
                r=random.random()
                if r < up_prob: 
                    row_vec.append(1)
                else:
                    row_vec.append(-1)
        M2d.append(row_vec)
    return np.array(M2d)           

def make_spin_ncr(L, n_down, updown):
    M2d=np.ones([L, L])
    L2 = L**2
    n=0
    if n_down <= int(L2/2):
        pass
    else:
        n_down = L2 - n_down
    while n < n_down:
        i, j = np.random.randint(L, size=2)
        if M2d[i][j] == 1:
            M2d[i][j] = -1
            n+=1
        else:
            pass

    if updown == 'up':
        return np.array(M2d)
    elif updown == 'down':
        M2d *= -1
        return np.array(M2d)


def periodicize(magnets):
    n=len(magnets)          # returns the first rank of the matrix
    m=len(magnets[0])
    if n != m:
        print("lattice is not square: %d by %d" % (n, m))
        exit(1)
    ss=np.zeros((n+2,n+2))
    for i in xrange(0,n):
        for j in xrange(0,n):
            ss[i+1][j+1]=magnets[i][j]
        ss[i+1][0]  = ss[i+1][n]
        ss[i+1][n+1]= ss[i+1][1]
    for i in xrange(0,n):
        ss[0][i+1]  = ss[n][i+1]
        ss[n+1][i+1]= ss[1][i+1]
    return ss

def total_magnetization(magnets_p):
    n=len(magnets_p)-2
    m=0
    for i in xrange(1,n+1):
        for j in xrange(1,n+1):
            m+=magnets_p[i][j]
    return m

def energy_cell_p(i,j,magnets_p):
    ii=i+1
    jj=j+1
    mag_sum = magnets_p[ii+1][jj]+magnets_p[ii-1][jj]+magnets_p[ii][jj+1]+magnets_p[ii][jj-1]
    ene = (-0.5)*magnets_p[ii][jj]*mag_sum
    return ene

def total_energy(magnets_p):
    n=len(magnets_p)-2
    e=0.
    for i in xrange(0,n):
        for j in xrange(0,n):
            ene = energy_cell_p(i,j,magnets_p)
            e += ene
    return 0.5*e    # total energy if half of sum of all the cell energies

def gen_data(N, L, sampling, outfile):

    #ofilename = ofile + ".npy"
    spin_length = L**2
    print "size of edge: ", L
    print "output data: ", outfile
    print "output formatt: spins, total energy"
    tf_data = []
    nsamples=[]
    if sampling == 'random':
        Nsample = N
        print "number of data: ", N
        for i in range(Nsample):
            mag_2d  = make_spin(L, sampling)
            mag_2dp = periodicize(mag_2d) 
            e=total_energy(mag_2dp)
            m=total_magnetization(mag_2dp)
            mag_1d = mag_2d.reshape(spin_length)
            dat = np.append(mag_1d,e)   # np append
            tf_data.append(dat)         # python append
    elif sampling == 'test_set' or sampling == 'nlog_scan':
        for i in range(0, spin_length+1):
            if sampling == 'test_set':
                pre = 1
            elif sampling == 'nlog_scan':
                pre = i
            n_C_r  = ncr(spin_length, i)
            ncr_disc = pre * np.log(float(n_C_r))
            num      = max(10, int(ncr_disc))
            nsamples.append(num)

        print "size of nsamples, number of data: ", len(nsamples), np.sum(nsamples)
        i = -1
        for num in nsamples:
            i += 1                  # start from i=0, 
            for n in range(num):
                if i <= 50:
                    mag_2d  = make_spin_ncr(L, i, 'up')   # for random selection of site (i,j)
                else:
                    mag_2d  = make_spin_ncr(L, i, 'down')   # for random selection of site (i,j)
                    
                mag_2dp = periodicize(mag_2d) 
                e=total_energy(mag_2dp)
                m=total_magnetization(mag_2dp)
                mag_1d = mag_2d.reshape(spin_length)
                ### check for spin in lattice
                '''
                if i == 100 :
                    print mag_2d
                    print e, m
                    exit(10)
                '''
                #print mag_1d
                dat = np.append(mag_1d,e)   # np append
                #print dat                   
                tf_data.append(dat)         # python append
                if i==0 or i==len(nsamples)-1:
                    break
    else:
        print "No sampling method"
        exit(1)
    
    np.save(outfile, tf_data)
    return 0

def main():
    parser = argparse.ArgumentParser(description='making spin configuration and add total energy')
    parser.add_argument('-L','-l','--length', default=10, type=int, help='side length of square of lattice')
    parser.add_argument('-N', '-n', '--Ndata', default=100000, type=int, help='number of data')
    parser.add_argument('-s', '-S', '--sample', choices={'nlog_scan', 'test_set', 'random'}, help='spin configuration')
    parser.add_argument('-o', '-O', '--outf', help='output filename')
    args = parser.parse_args()

    filename = 'spin2dL'+str(args.length)
    if not args.sample:
        filename += '_try' 
    else :
        filename += '_' + args.sample
    #filename += '.dat'
    if not args.outf:
        pass
    else:
        filename = args.outf

    gen_data(args.Ndata, args.length, args.sample, filename)

if __name__ == "__main__":
    main()
